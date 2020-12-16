.PHONY: make_conda update_conda remove_conda update_data update_data_dry

.ONESHELL:

PROJECT_NAME = nmr_assign_state 

make_conda:
	conda env create -f environment.yml
	ipython kernel install --user --name=$(PROJECT_NAME)

update_conda:
	conda env update --file environment.yml
	ipython kernel install --user --name=$(PROJECT_NAME)

remove_conda:
	conda remove --name=$(PROJECT_NAME) --all

update_data:
	 while read a b;do rsync -rauLih --progress  --include-from=data/raw/include_list.txt   $$a data/raw/$$b;done < data/raw/dir_list.txt
	 python  src/data/get_external_data.py

update_data_dry:
	 while read a b;do rsync -raunLih  --include-from=data/raw/include_list.txt   $$a data/raw/$$b;done < data/raw/dir_list.txt

help:
	@echo "Possible options:"
	@echo "make_conda"
	@echo "update_conda"
	@echo "remove_conda"
	@echo "update_data"
	@echo "update_data_dry"
