def get_chemical_shifts(univ,temp_dir='./',split_size=100,skip=1,method='sparta_plus'):
    import MDAnalysis as mda
    import mdtraj as md
    import numpy as np
    import pandas as pd
    import os
    from tqdm import tqdm
    import warnings
    warnings.filterwarnings("ignore")

    '''
    Use sparta_plus with the interface of mdtraj to get the chemical shifts of     the trajectory. The segids are concatenated as independent trajectories.

    Parameters
    ----------
    univ : Universe object that contains the protein.
    temp_dir : Temporary directory for the tmp trajectory.
    split_size : How to split the trajectory so the /tmp does not fill.
    skip : Amount of frames skiped.
      
    Returns
    -------
    df: Dataframe of chemical shifts (n_cs_nuclei,n_frames*n_segids)
    '''
    assert isinstance(temp_dir,str), f'{temp_dir} should be a str'
    assert isinstance(split_size,int), f'{split_size} should be a int'
    assert isinstance(skip,int), f'{skip} should be a int'
    
    sel_protein=univ.select_atoms('protein')
    xtc = f"{temp_dir}/tmp.xtc"
    pdb = f"{temp_dir}/tmp.pdb"
    
    
    #Which method:
    if method == 'sparta_plus': 
        method_function = md.nmr.chemical_shifts_spartaplus
    elif method == 'ppm': 
        method_function = md.nmr.chemical_shifts_ppm
    else:
        raise Exception(f'Method {method} not recognized.')
    
    list_of_splits=[]
    for segid in np.unique(sel_protein.segids):
        sel_segid=sel_protein.select_atoms(
            f'segid {segid}') 
        with mda.Writer(xtc, n_atoms=sel_segid.atoms.n_atoms) as W:
            for ts in univ.trajectory[::skip]:
                W.write(sel_segid)
        sel_segid.write(pdb)
        
        #Load data and calculate NMR-CS
        trj= md.load(xtc,top=pdb)
        n_frames=trj.n_frames
        for j in tqdm(range(0,n_frames,split_size),
                      desc=f't({segid})',
                      leave=True,smoothing=1):
            list_of_splits.append(
           method_function(trj[j:j+split_size])
            )
        df = pd.concat(list_of_splits,axis=1)
        df.columns=np.arange(0,len(df.columns))       
    os.remove(xtc)
    os.remove(pdb)   
    return df
