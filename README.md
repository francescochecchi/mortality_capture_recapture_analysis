## Generic R scripts for capture-recapture analysis of multiple lists of deceased persons
## Explanatory note

October 2021

Author:
Francesco Checchi (London School of Hygiene and Tropical Medicine)

francesco.checchi@lshtm.ac.uk

Contributors: Rahaf Abu Koura, Mervat Alhaffar, Maysoon Dahab, Pasquale Pepe

Donor support: United Kingdom Foreign and Commonwealth Office, British Council


## Background
Capture-recapture analysis (also known as multiple systems estimation or mark-recapture) is a well-established statistical technique to infer the true number of 'cases' from two or more lists, each of which may be incomplete to an unknown extent. The technique was developed for ecological studies but has been applied in multiple fields, including epidemiology. A key data requirement for capture-recapture analysis is that each list contain enough identifier information on the 'case' to enable confident record linkage, i.e. determination of whether any case also appears on other lists.

This set of R scripts was written in the context of studies of population mortality in settings with incomplete vital registration (Sudan, Yemen). It allows for analysis of two, three or four lists of individual decedents, composed by asking key community informants to recall deaths within their community during a specified retrospective period. Apart from the ubiquitous reference to deaths, however, the scripts could be used to analyse lists of any other event (e.g. cases of a given disease, persons who experienced a specific kind of injury, etc.).

## Description of the R scripts
* script `control_code` installs required R packages and identifies the working directory; reads in datasets, as well as user-defined parameters governing how analysis should be done for each location; cleans each dataset and prepares any strata of interest (age, gender, period); compares the lists in terms of key characteristics, and visualises their overlap; and performs capture-recapture analysis (see below). The script generates reasonably well-formatted tables and figures, which will appear on the same working directory;
* script `bespoke_functions` contains a number of user-defined functions for performing the above cleaning, descriptive and analytic steps.
Unfortunately the scripts do not yet perform record linkage. This will have to be done by the user ahead of analysis, using manual and/or automated techniques.

## How to implement analysis
The user should prepare an Excel file whose name should include the string `list_data` (e.g. `list_data_studya.xlsx`), structured as per the dummy file included in this directory. The *dictionary* worksheet (tab) of the Excel file contains more specifications about each data column. Critically, the dataset must be organised so that one line is one unique decedent, who may appear on only one list or several. Furthermore, the dataset from each study (or location) should be included as its own worksheet, with the worksheet name corresponding exactly to the location variable in the <parameters> worksheet. The user interacts with the analysis through the *parameters* worksheet, by specifying all the meta-variables and parameters listed in the *dictionary*. These allow for specifying, for each site / study, a variable number of lists (between two and four), plausibility criteria for candidate models (see below), restrictions to the dataset in terms of gender, age or period, cut-offs to construct strata and/or categorical variables, an exposure variable and variables to adjust for. Key variables (e.g. dates) should be entered using the column names and format specified in the *directory* worksheet; there is a bit of flexibility built in, e.g. the script will recognise both numeric and text versions of months and different special characters to separate day, month and year. The script will also automatically define eligibility of each observation, based on availability of minimal data. On the other hand, if the user specifies ‘impossible’ sets  of parameters (e.g. the age range is restricted to 15y and above but age category cut-offs include the under 5y), the script will either throw out a warning on the console, or stop executing, with an error message.

The Excel file should be placed in the same directory as the two R scripts. The `control_code` script will automatically recognise the directory. The user needn't alter any of the scripts. To run analysis, the scripts should be executed from R itself (freely downloadable from https://www.r-project.org/) or (preferably) from RStudio (https://www.rstudio.com/products/rstudio/download/: the free desktop version suffices, but note that RStudio also requires R to be installed). Simply open the `control_code.r` script and press Ctrl-Alt-R to run the analysis. All output graphs and tables should appear automatically on the directory. For each location, analysis should take no more than a few minutes on a standard laptop.

## Details on the capture-recapture estimates
If the dataset only contains two lists, a simple analysis will be done using the Chapman estimator and associated 95% confidence interval (see for ex. https://en.wikipedia.org/wiki/Mark_and_recapture). Two-list capture-recapture analysis has several limitations, specifically that the two lists must be independent of each other. This assumption can be relaxed if more than two lists are collected, by fitting a variety of candidate log-linear Poisson models that include various potential interaction terms representing the (mutual) dependence of one list on another (for two-way interactions) or on other lists (for higher-order interactions). Two-list analysis will ignore any exposure or adjustment variables that may have been specified, but analysis for each desired stratum will be performed.

The scripts allow for three- or four-list analysis (beyond this, the number of candidate models becomes unmanageably large). In both instances, each candidate log-linear model (featuring different combinations of interaction terms among lists, and any exposure/confounder variables specified) is fit to an expanded dataset featuring each possible list outcome for each decedent, and values of the m000 or m0000 cells (i.e. deaths that are not on any of the lists) are predicted; models that produce implausible m000(0) estimates are screened out based on a user-defined plausibility parameter (i.e. the maximum tolerated ratio of the total estimated to the total observed deaths), and a cut-off for the model’s likelihood ratio p-value (compared to the most parsimonious model) is also applied, above which models are excluded on suspicion of overfitting the data. For the parameterisation of the analysis, see Rossi et al. (https://rivista-statistica.unibo.it/article/view/9854).

Instead of selecting a single candidate model, the script then follows Rossi et al.’s (https://rivista-statistica.unibo.it/article/view/3593/2945) method for model averaging, whereby a mean estimate of m000(0) is arrived at by weighting each short-listed model based on its Akaike Information Criterion, itself related to the model’s goodness-of-fit and relative parsimony. The sensitivity of each list and all lists combined are also computed.

## Notes of caution
The following should be borne in mind when selecting parameters for the analysis, and interpreting output:
* The sparseness of the data should be considered carefully, for example by running the descriptive part of the analysis with a given set of stratification cut-offs, and observing the extent of overlap among lists: generally, situations in which several cells in the k-way contingency table are empty ( = 0) will result in poor or no model fit. The output will show which, if any, of the candidate models actually fit successfully. While the analysis will complete even if few (or only one) of the candidate models are short-listed to the model averaging step, the user should look at the output carefully, and avoid relying on estimates that are the result of only a few candidate models out of all the possible ones. Unless more data can be collected, the only mitigation for sparseness at this stage is to minimise the analysis strata (e.g. age cut-offs), and/or restrict analysis to a period, age range or gender for which there appears to be sufficiently rich data.
* Relatedly, four-list analysis, while theoretically preferable to three lists, will generally be better supported by fairly rich data, e.g. reasonably large datasets with cases falling within each of the contingency table cells. If data are sparse, one strategy could be to reduce the system to three lists only, either by excluding one of the lists or by combining two lists into one (if and only if this is reasonable in practice, e.g. it might be OK to combine reports of cases from informal pharmacies with those from formal health facilities, but religious leaders should probably not be combined with formal health facilities, as they would usually be very different sources).
* The output of each candidate model shows results when each or all of the adjustment variables (confounders) are removed from the model formula: generally, confounders should only be retained if they appreciably modify the m000(0) estimates. At least with sparse data, models featuring confounders may also fit more poorly or not at all, another reason to introduce confounders with caution and inspect output carefully.
