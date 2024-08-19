PFAS Chemical Categorisation
==============================

Code repository supporting manuscript titled "Development of Chemical Categories for Per- and Polyfluoroalkyl Substances (PFAS) and the Proof-of-Concept Approach to the Identification of Potential Candidates for Tiered Toxicological Testing and Human Health Assessment". 


Please cite: Patlewicz G., Judson R., Williams A. J., Butler T., Barone Jr. S ., Carstens K. E., Cowden J., Dawson J. L., Degitz S. , Fay K., Henry T. R., Lowit A.,  Padilla S., Paul Friedman K., Phillips M. B., Turk D., Wambaugh J., Wetmore B., Thomas R.S. Development of Chemical Categories for Per- and Polyfluoroalkyl Substances (PFAS) and the Proof-of-Concept Approach to the Identification of Potential Candidates for Tiered Toxicological Testing and Human Health Assessment. *Computational Toxicology* **2024** https://doi.org/10.1016/j.comtox.2024.100327

Supplementary information referenced in the manuscript are available as a compressed tar file.
Raw data files underpinning the analysis are available as a compressed tar file. 

These files can be accessed at doi.org/10.23645/epacomptox.26524327

Code repository is provided on an "as is" basis. 

Project Organization
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── docs               <- A default Sphinx project; see sphinx-doc.org for details
    │
    ├── models             <- Trained and serialized models, model predictions, or model summaries
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so src can be imported
    ├── src                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── data           <- Scripts to download or generate data
    │   │   └── make_dataset.py
    │   │
    │   ├── features       <- Scripts to turn raw data into features for modeling
    │   │   └── build_features.py
    │   │
    │   ├── models         <- Scripts to train models and then use trained models to make
    │   │   │                 predictions
    │   │   ├── predict_model.py
    │   │   └── train_model.py
    │   │
    │   └── visualization  <- Scripts to create exploratory and results oriented visualizations
    │       └── visualize.py
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
