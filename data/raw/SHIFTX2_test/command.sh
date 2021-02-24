conda activate nmr_assign_state_for_shiftx2
shiftx2.py -i protein_3FB5_lb.pdb -c A -e -p 4 -t 290 -v
sort -k1 --field-separator="," -g  protein_3FB5_lb.pdb.cs|grep "\(,CA,\|,CB,\|,C,\|,N,\)" > protein_3FB5_lb.pdb.cs_sorted
sort -k1 --field-separator="," -g  protein_3FB5_lb.pdb.cf.sxp|grep "\(,CA,\|,CB,\|,C,\|,N,\)" > protein_3FB5_lb.pdb.cf.sxp_sorted
