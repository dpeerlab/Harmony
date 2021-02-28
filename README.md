Harmony
------

Harmony is a unified framework for data visualization, analysis and interpretation of scRNA-seq data measured across discrete time points.Harmony constructs an augmented affinity matrix by augmenting the kNN graph affinity matrix with mutually nearest neighbors between successive time points. This augmented affinity matrix forms the basis for generated a force directed layout for visualization and also serves as input for computing the diffusion operator which can be used for trajectory detection using Palantir


#### Installation and dependencies
1. Harmony has been implemented in Python3 and can be installed using:

        $> pip install harmonyTS
        $> pip install palantir

2. Harmony depends on a number of `python3` packages available on pypi and these dependencies are listed in `setup.py`
All the dependencies will be automatically installed using the above commands

3. To uninstall:
		
		$> pip uninstall harmonyTS

4. If you would like to determine gene expression trends, please install <a href="https://cran.r-project.org"> R <a> programming language and the R package <a href="https://cran.r-project.org/web/packages/gam/">GAM </a>. You will also need to install the rpy2 module using 
	
		$> pip install rpy2
		

#### Usage

A tutorial on Harmony usage and results visualization for single cell RNA-seq data can be found in this notebook: http://nbviewer.jupyter.org/github/dpeerlab/Harmony/blob/master/notebooks/Harmony_sample_notebook.ipynb

The datasets generated as part of the manuscript and harmozined using Harmony are available for exploration at: [endoderm-explorer.com](https://endoderm-explorer.com)


#### Citations
Harmony was used to harmonize datasets across multiple time points in our manuscript characterizing mouse gut endoderm development.  This manuscript is available at [Nature](https://www.nature.com/articles/s41586-019-1127-1). If you use Harmony for your work, please cite our paper.

