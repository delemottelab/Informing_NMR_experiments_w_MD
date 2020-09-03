.PHONY: make_conda update_conda remove_conda

.ONESHELL:

PROJECT_NAME = nmr_assign_state 

make_conda:
	conda env create -f environment.yml
	ipython kernel install --user --name=$(PROJECT_NAME)

update_conda:
	conda env update --file environment.yml

remove_conda:
	conda remove --name $(PROJECT_NAME) --all

