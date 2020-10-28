def get_chemical_shifts(univ,temp_dir='./',split_size=100,skip=1,method='sparta_plus',temperature=298.,pH=5.0, protein_selection = 'protein'):
    import MDAnalysis as mda
    import mdtraj as md
    import numpy as np
    import pandas as pd
    import os
    from tqdm import tqdm
    import warnings
    warnings.filterwarnings("ignore")

    '''
    Use sparta_plus with the interface of mdtraj to get the chemical shifts of the trajectory. The segids are concatenated as independent trajectories.
    The protein is assumed to have 4 identical subunits. 

    Parameters
    ----------
    univ : Universe object that contains the protein.
    temp_dir : Temporary directory for the tmp trajectory.
    split_size : How to split the trajectory so the /tmp does not fill.
    skip : Amount of frames skiped.
    pH: pH for shiftx2 calculations.
    temperature: temperature for shiftx2 calculations.
    protein_selection: protein selection for MDAnalysis.
      
    Returns
    -------
    df: Dataframe of chemical shifts (n_cs_nuclei,n_frames*n_segids)
    '''
    assert isinstance(temp_dir,str), {temp_dir} + ' should be a str'
    assert isinstance(protein_selection,str), 'protein_selection should be a str'

    assert isinstance(split_size,int), split_size +  ' should be a int'
    assert isinstance(skip,int), skip + ' should be a int'
    assert isinstance(temperature,float), str(temperature) + ' should be a float'
    assert isinstance(pH,float), str(pH) + ' should be a float'
    
    sel_protein=univ.select_atoms(protein_selection)
    resids = sel_protein.residues.resids.copy()
    gap = 10
    last_resid = -gap
    for i,seg in enumerate(univ.segments):
        seg.residues.resids = seg.residues.resids + gap + last_resid
        last_resid = seg.residues.resids[-1]
    renumber = sel_protein.residues.resids
    resids_dic = { i:j  for i,j in zip(renumber,resids)}
    sel_protein.segments.segids = 'A'
    
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
    
    
     
    with mda.Writer(xtc, n_atoms=sel_protein.atoms.n_atoms) as W:
        for ts in univ.trajectory[::skip]:
            W.write(sel_protein)
    sel_protein.write(pdb)

    #Load data and calculate NMR-CS
    trj= md.load(xtc,top=pdb)
    n_frames=trj.n_frames
    list_of_splits = []
    for j in tqdm(range(0,n_frames,split_size),
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
    index = pd.MultiIndex.from_arrays([[ resids_dic[i[0]] for i in df.index ], [ i[1] for i in df.index ]], names = ['resid', 'nuclei'])
    df.index = index
    df = df.T
    df = df.loc[:,pd.IndexSlice[:,['N', 'CA', 'CB', 'C']]]
    dictionary = {}
    for col in np.unique(df.columns):
        k = df.loc[:,col]
        if k.shape[1] == 4:
            dictionary[col] = np.array(k).reshape(k.shape[0]*k.shape[1])
    df = pd.DataFrame(dictionary)
    df = df.rename_axis(['resid', 'nuclei'], axis=1)
    os.remove(xtc)
    os.remove(pdb)
    return df
