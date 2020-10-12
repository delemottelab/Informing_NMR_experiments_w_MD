def get_chemical_shifts(univ,temp_dir='./',split_size=100,skip=1,method='sparta_plus',temperature=298.,pH=5.0):
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
    pH: pH for shiftx2 calculations.
    temperature: temperature for shiftx2 calculations.
      
    Returns
    -------
    df: Dataframe of chemical shifts (n_cs_nuclei,n_frames*n_segids)
    '''
    assert isinstance(temp_dir,str), {temp_dir} + ' should be a str'
    assert isinstance(split_size,int), split_size +  ' should be a int'
    assert isinstance(skip,int), skip + ' should be a int'
    assert isinstance(temperature,float), str(temperature) + ' should be a float'
    assert isinstance(pH,float), str(pH) + ' should be a float'
    
    sel_protein=univ.select_atoms('protein and not resid 26 121')
    xtc = temp_dir+"/tmp.xtc"
    pdb = temp_dir+"/tmp.pdb"
    
    
    #Which method:
    if method == 'sparta_plus': 
        method_function = md.nmr.chemical_shifts_spartaplus
    elif method == 'ppm': 
        method_function = md.nmr.chemical_shifts_ppm
    elif method == 'shiftx2': 
        method_function = md.nmr.chemical_shifts_shiftx2
    else:
        raise Exception('Method' + method + 'not recognized.')
    
    list_of_splits=[]
    for segid in np.unique(sel_protein.segids):
        sel_segid=sel_protein.select_atoms(
            'segid '+segid) 
        with mda.Writer(xtc, n_atoms=sel_segid.atoms.n_atoms) as W:
            for ts in univ.trajectory[::skip]:
                W.write(sel_segid)
        sel_segid.write(pdb)
        
        #Load data and calculate NMR-CS
        trj= md.load(xtc,top=pdb)
        n_frames=trj.n_frames
        for j in tqdm(range(0,n_frames,split_size),
                      desc='t('+segid+')',
                      leave=True,smoothing=1):
            if method == 'shiftx2':
                list_of_splits.append(
               method_function(trj[j:j+split_size],pH=pH,temperature=temperature)
                )
            else:
                list_of_splits.append(
               method_function(trj[j:j+split_size])
                )
        df = pd.concat(list_of_splits,axis=1)
        df.columns=np.arange(0,len(df.columns))       
    os.remove(xtc)
    os.remove(pdb)   
    return df