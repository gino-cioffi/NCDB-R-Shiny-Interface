# NCDB-R-Shiny-Interface

This application was created to show basic descriptive statistics and overall survival analyses using National Cancer Database (NCDB) Data.

### Tab Panel One: Descriptive Statistics

Panel one allows the user to visually see the distribution of selected variables among the selected population. 

  If the user selects a categorical variable (e.g Sex) they will be presented with a bar chart. The user can select an additional categorical variable as a sub-plot, which will show the distribution of this variable within each level of the initially selected variable and a side-by-side bar chart. The user can select to display proportions instead of counts on the bar chart as well.
  A table of counts and percentages is generated below the figure. If the user selects a subplot, a Chi-Square test is performed between both selected variables and the associated P value is presented in the last column (it is assumed that the dataset is large enough that all cell counts > 5, if not P value may not be interpretable). 
  
  If the user selects a continuous variable (e.g. Age) they will be presented with a histogram. The user can select an additional categorical variable as a sub-plot, which will show the distribution of the continuous variable within each level of the categorical variable as a separate histogram. 
  A table of means and standard deviations is generated below the figure. If the user selects a subplot, an analysis of variance (ANOVA) is performed to compared mean values across levels and the associated P value is presented in the last column (similarly, due to the size of the dataset it is assumed that the distributions of the data meet the assumption of the statistical test, if not the P value may not be interpretable). 
  
### Tab Panel Two: Survival Analysis - Kaplan-Meier Survival Graph

Panel two allows the user to visually see survival differences between levels of selected variables among the selected population. 

  The user can select a categorical variable to display a Kaplan-Meier survival curve. The user can use the slider bar to select, to see survival probability for each level at any given point in time, which appears in a table beneath the figure.  Median survival and the 95% confidence interval for each level is always present at the bottom table on the page.
