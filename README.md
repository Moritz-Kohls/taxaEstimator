### taxaEstimator
Estimating the taxonomic composition of viral sequences in a biological sample processed by next-generation sequencing is an important step for comparative metagenomics. For that purpose, sequencing reads are usually classified by mapping them against a database of known viral reference genomes. This fails, however, to classify reads from novel viruses and quasispecies whose reference sequences are not yet available in public databases.

In order to circumvent these problems, the feasibility and performance of an artificial neural network to classify sequencing reads to taxonomic orders is studied. For that purpose, taxonomy and genomic data from the NCBI database are used to sample artificial reads from viruses with known taxonomic attribution. From these reads, input features are derived by a feature selection method. Based on the computed training data, an artificial neural network is fitted and applied to classify single viral read sequences to different taxa. While we observe many misclassifications with this neural network, we employ additionally simulated test samples with known taxonomic memberships to correct the estimation of taxa distributions in the final test sample. Therefore, new formulas are introduced to statistically estimate taxa frequencies.

Prediction accuracy of the fitted models is evaluated on the simulated test data and classification results are summarised in a confusion matrix, from which sensitivity, specificity as well as predictive values are calculated. The prediction accuracy of the artificial neural network is considerably higher than for random classification and posterior (i.e. corrected) estimation of taxa frequencies is closer to the true distribution in the test data than simple classification or mapping results.

The general principle of using additionally simulated test samples with known class memberships could also be helpful to correct class distributions obtained by other machine learning applications.

### Installation
If you want to install R-packages from GitHub, at first you need to install the R-package **devtools** from CRAN:

```r
install.packages("devtools")
library(devtools)
```
After this, you can install packages from GitHub by using the R-function ```install_github``` and specifying the author and package name. For example, if you want to install the R-package **taxaEstimator**, just type the following:

```r
install_github("Moritz-Kohls/taxaEstimator")
```

Type `help(package = taxaEstimator)` to open the package's help page.
