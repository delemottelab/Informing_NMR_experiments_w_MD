'''
Some functions to handle xr data.
'''
def combine_posterior_variables(posterior, variables_to_combine, new_name_var,
                                new_dim_name, new_dim_names):
    import xarray as xr
    '''
    Combine several variables of a posterior to a single one
    with new dimensions.

    Parameters
    ----------
    posterior: xarray posterior.
    variables_to_combine: list of variables names.
    new_name_var: str, name of the new combined variable.
    new_dim_name: New dimension name.
    new_dim_names: list of new_dim names strings
    Returns
    -------
    Posterior with combined dimension.
    '''
    posteriors = (posterior[variable] for variable in variables_to_combine)
    posterior[new_name_var] = xr.concat(
        posteriors, dim=new_dim_name).assign_coords(
        {new_dim_name: new_dim_names})
    posterior = posterior.drop(variables_to_combine)
    return posterior