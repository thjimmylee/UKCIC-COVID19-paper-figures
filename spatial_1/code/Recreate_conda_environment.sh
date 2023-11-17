conda env create -f ../environment.yml

conda activate spacejam

pip install --no-cache-dir "pymc3>=3.8,<3.10"

pip install --no-cache-dir git+https://github.com/BayraktarLab/cell2location.git

pip install --no-cache-dir git+https://github.com/BayraktarLab/CountCorrect/

python -m ipykernel install --user --name spacejam --display-name "Python (SpaceJam)"