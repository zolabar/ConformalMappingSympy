# ConformalMappingSymPy 

This repository hosts code that illustrates the content of the presentation and proceeding *Conformal Mappings with Sympy: Towards Python-driven Analytical Modeling in Physics* a collaborative work of [Lauer-Baré Z.](https://orcid.org/0000-0002-7083-6909) and [Gaertig E.](https://orcid.org/0000-0003-1728-6466) presented on THE 20th PYTHON IN SCIENCE CONF. (SCIPY 2021).

Python scripts and jupyter notebooks are provided.

Please refer to 

Lauer-Baré Z.and Gaertig E., *Conformal Mappings with SymPy: Towards Python-driven Analytical Modeling in Physics*, PROC. OF THE 20th PYTHON IN SCIENCE CONF. (SCIPY 2021)

when using formulae or code from this repository.

The theoretical methods used here are conformal mappings, inspired by [PHW33](https://www.tandfonline.com/doi/abs/10.1080/14786443309462212) and [BC09](https://www.mheducation.com/highered/product/complex-variables-applications-brown-churchill/M9780073383170.html) and Taylor-expansions, following [LGK21](https://journals.riverpublishers.com/index.php/IJFP/article/view/5535).

## Transformation of eccentric annulus to concentric annulus

```code_block_moebius.py``` and ```moebius.ipynb``` with a Möbius transform of the type

<img src="https://render.githubusercontent.com/render/math?math=w(z)=\frac{z %2B ia}{az %2B i}">

## Transformation of eccentric annulus to rectangle

```code_block_bipolar.py``` and ```bipolar.ipynb``` with a conformal mapping related to bipolar coordinates

<img src="https://render.githubusercontent.com/render/math?math=w(z)=2\cdot \tan^{-1}\left(\frac{z %2B i\gamma}{c}\right)">

## Postprocessing

The postprocessing is shown in the file ```moebius.ipynb``` , due to LaTeX rendering of web browser based *jupyter notebook*.

Flow force calculation with ```diff``` and Taylor expansion of force in the gap 

<img src="https://render.githubusercontent.com/render/math?math=\delta">

with ```series```.

## Literature

[BC09] [Brown JW, Churchill RV. Complex variables and applications eighth edition, McGraw-Hill Book Company; 2009](https://www.mheducation.com/highered/product/complex-variables-applications-brown-churchill/M9780073383170.html)

[LGK21] [Lauer-Baré Z., Gaertig E., Krebs J., Arndt C., Sleziona A., Gensel A. A note on leakage jet forces: Application in the modelling of digital twins of hydraulic valves, International Journal of Fluid Power, 2021, Vol. 22 (1), 113–146](https://journals.riverpublishers.com/index.php/IJFP/article/view/5535)

[PHW33] [N.A.V. Piercy D.Sc., M.S. Hooper & H.F. Winny Ph.D. LIII. Viscous flow through pipes with cores, The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science, 1933](https://www.tandfonline.com/doi/abs/10.1080/14786443309462212)
