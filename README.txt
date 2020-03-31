This is the version which gets deployed online. Please make edits in the development version "shiny_dev".

To load it locally:
* make sure you've got the R package "shiny". It may also force you to install some other packages but hopefully the errors will be informative
* load R
* set your working directory to here
* do "shiny::runApp()"

Parameters start off at the values saved in "initial_params.RData", which was generated from "initial_params.R", but these are now hard-coded in to match the paper.