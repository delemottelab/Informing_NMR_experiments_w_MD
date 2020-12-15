import matplotlib.pyplot as plt
import arviz as az
RANDOM_SEED = 281090

raw_data_dir = '../data/raw/'
interim_data_dir = '../data/interim/'
processed_data_dir = '../data/processed/'
external_data_dir = '../data/external/'
models_dir = '../models/bayesian/'
reports_dir = '../reports/'
model_name = 'skew_model'

method_to_name = {
    'sparta_plus': 'SPARTA+',
    'shiftx2': 'SHIFTX2'
}
nuclei_to_name = {
    'CA': r'C$_\alpha$',
    'CB': r'C$_\beta$',
    'C': 'C',
    'N': 'N',
}
state_to_name = {
    'like_o' : 'Open',
    'like_fo' : 'Fully Open',
    'like_c' : 'Closed',
}
for nucleus in nuclei_to_name.keys():
    for method in method_to_name.keys():
        for state in state_to_name.keys():
            print(nucleus,method,state)
            model_path = models_dir + f"{model_name}_{method}_{nucleus}.nc"
            my_model = az.from_netcdf(model_path)
            resids = my_model.posterior.resid
            n = resids.shape[0]
            fig, ax = plt.subplots(n // 6 + 1, 6, figsize=(13,15))
            for i in range(6 - n % 6):
                fig.delaxes(ax[-1,-i-1])
            ax = fig.axes
            with az.rc_context(rc={'plot.max_subplots': None}):
                az.plot_ppc(my_model, flatten=['step'], var_names = [state], random_seed=RANDOM_SEED, ax=ax)
            for r, a in zip(resids.to_index(), ax):
                a.set_title(f'{r}', size=12)
                a.set_xlabel('')
                a.legend_.set_visible(False)
            fig.suptitle(f'Posterior Predictive Check {method_to_name[method]} {nuclei_to_name[nucleus]} {state_to_name[state]} (ppm)', y =1.0, size =20)
            fig.tight_layout()
            plt.savefig(f'{reports_dir}for_print/ppc_{state_to_name[state]}_{nucleus}_{method}.png')
