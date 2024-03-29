[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zolabar/ConformalMappingSympy/HEAD)
# Conformal Mappings with SymPy 

This repository hosts code that illustrates the content of the presentation and conference proceeding *Conformal Mappings with SymPy: Towards Python-driven Analytical Modeling in Physics*, a collaborative work of [Zoufiné Lauer-Baré](https://orcid.org/0000-0002-7083-6909) and [Erich Gaertig](https://orcid.org/0000-0003-1728-6466) presented on [THE 20th PYTHON IN SCIENCE CONF. (SCIPY 2021)](https://www.scipy2021.scipy.org/).

Python scripts and Jupyter notebooks are provided.

Please refer to 

Lauer-Baré Z. and Gaertig E., [*Conformal Mappings with SymPy: Towards Python-driven Analytical Modeling in Physics*. Lauer-Baré, Z. & Gaertig, E. In Agarwal, M., Calloway, C., Niederhut, D., & Shupe, D., editors, Proceedings of the 20th Python in Science Conference, pages 85 - 93, 2021. ](https://conference.scipy.org/proceedings/scipy2021/lauer_bare_gaertig.html)

when using formulae, code, figures or animations from this repository. The conference talk can be seen on the [Enthought YouTube channel](https://www.youtube.com/watch?v=P5ybpjv2uDA).

The theoretical methods used here are conformal mappings, inspired by [PHW33](https://www.tandfonline.com/doi/abs/10.1080/14786443309462212) and [BC09](https://www.mheducation.com/highered/product/complex-variables-applications-brown-churchill/M9780073383170.html) and Taylor-expansions, following [LGK21](https://journals.riverpublishers.com/index.php/IJFP/article/view/5535) and [LGKS23](https://publications.rwth-aachen.de/record/957191). These methods are used to solve the Stokes problem in an eccentric annular domain for Couette-Poiseuille flow and to calculate the corresponding flow force in a postprocessing step, as well as analyzing the limits for small gaps. The context of this work is the modelling of viscous fluid power systems (see [LGK21](https://journals.riverpublishers.com/index.php/IJFP/article/view/5535) for more details).

Applications of conformal mappings with SymPy in the context of inviscid irrotational flow can be found on [Plotting streamlines with Matplotlib and SymPy](https://tonysyu.github.io/plotting-streamlines-with-matplotlib-and-sympy.html#.YPf_rKjwhPb) (T. S. Yu).
Further applications of conformal mappings with SymPy in the context of inviscid irrotational flow applied to naval engineering are described in [G21](https://www.mdpi.com/2077-1312/9/2/108) with an open Python code repository in [naval Python and SymPy](https://zenodo.org/record/4452633#.YPpnYegzZPZ).

Further, a Python package with cloud computing possibility via binder, for visualizing conformal mappings interactively, based on SymPy, NumPy and Plotly can be found at the page  [conformalMaps](https://github.com/im-AMS/Conformal-Maps); see [LA21](https://doi.org/10.5281/zenodo.5717868).

## Transformation of eccentric annulus to concentric annulus

```code_block_moebius.py``` and ```moebius.ipynb``` with a Möbius transform of the type

<img src="https://render.githubusercontent.com/render/math?math=w(z)=\frac{z %2B ia}{az %2B i}">

The following animation was created with an adapted version of the code from [the interactive Python code for conformal mappings mentioned above](https://github.com/im-AMS/Conformal-Maps), where the results from [LG21](https://conference.scipy.org/proceedings/scipy2021/lauer_bare_gaertig.html) were implemented.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zolabar/Conformal-Maps/HEAD)
![](Figures/moebius.gif)

## Transformation of eccentric annulus to rectangle

```code_block_bipolar.py``` and ```bipolar.ipynb``` with a conformal mapping related to bipolar coordinates

<img src="https://render.githubusercontent.com/render/math?math=w(z)=2\cdot \tan^{-1}\left(\frac{z %2B i\gamma}{c}\right)">

The following animation was created with an adapted version of the code from [the interactive Python code for conformal mappings mentioned above](https://github.com/im-AMS/Conformal-Maps), where the results from [LG21](https://conference.scipy.org/proceedings/scipy2021/lauer_bare_gaertig.html) were implemented.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/zolabar/Conformal-Maps/HEAD)
![](Figures/mapping_arctan_colored_boundary.gif)

## Postprocessing

The postprocessing is shown in the file ```moebius.ipynb``` , due to LaTeX rendering of web browser based *Jupyter notebook*.

Flow force calculation with ```diff``` and Taylor expansion of force in the gap 

<img src="https://render.githubusercontent.com/render/math?math=\delta">

with ```series```.

## Literature

[BC09] [Brown J.W., Churchill R.V., Complex variables and applications, Eighth edition, McGraw-Hill Book Company; 2009](https://www.mheducation.com/highered/product/complex-variables-applications-brown-churchill/M9780073383170.html)

[G21] [Grm, A., Ships Added Mass Effect on a Flexible Mooring Dolphin in Berthing Manoeuvre. J. Mar. Sci. Eng. 2021, 9, 108. https://doi.org/10.3390/jmse9020108](https://www.mdpi.com/2077-1312/9/2/108)

[LGK21] [Lauer-Baré Z., Gaertig E., Krebs J., Arndt C., Sleziona C., Gensel A., A note on leakage jet forces: Application in the modelling of digital twins of hydraulic valves, International Journal of Fluid Power, 2021, Vol. 22 (1), 113–146](https://journals.riverpublishers.com/index.php/IJFP/article/view/5535)

[LGKS23] [Lauer-Baré Z., Gaertig E., Krebs J., Sleziona C. Analytical formulae in fluid power, quo vadis in times of CFD and I4.0? 13th International Fluid Power Conference, p. 866-879, 2023](https://publications.rwth-aachen.de/record/957191)

[LG21] Lauer-Baré Z. and Gaertig E., [*Conformal Mappings with SymPy: Towards Python-driven Analytical Modeling in Physics*. Lauer-Baré, Z. & Gaertig, E. In Agarwal, M., Calloway, C., Niederhut, D., & Shupe, D., editors, Proceedings of the 20th Python in Science Conference, pages 85 - 93, 2021. ](https://conference.scipy.org/proceedings/scipy2021/lauer_bare_gaertig.html)

[LA21] Lauer-Baré Z. and Aditya, [*Conformal-Maps: Code for interactive conformal mapping with python and jupyter notebook* (v1.0.1). Zenodo. 2021 https://doi.org/10.5281/zenodo.5717868.](https://doi.org/10.5281/zenodo.5717868)

[PHW33] [Piercy N.A.V., Hooper M.S., Winny H.F., LIII. Viscous flow through pipes with cores, The London, Edinburgh, and Dublin Philosophical Magazine and Journal of Science, 1933](https://www.tandfonline.com/doi/abs/10.1080/14786443309462212)
