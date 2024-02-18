# Selection_and_aggregation_of_conformal_prediction_sets


This repo contains the code and data for the simulation and real data application results in [this paper](https://arxiv.org/abs/2104.13871).  

- ```methods_functions.R``` contains all the functions used in this work. `efcp.fun` is the function to implement Efficiency First Conformal Prediction(EFCP), `vfcp.fun` is for implementing Validity First Conformal Prediction(VFCP), `conf_CQR_reg` is used to implement Conformalized Quantile Regression(CQR).

The files that have names ending with '_simulation.R' uses the functions in ```methods_functions.R``` and are used to generate the results (mostly summary statistics such as coverage and size of prediction intervals contained in ```.RData``` files). The files ending with ```_plot.R``` are used to yield the plots illustrated in the paper. 

##
For example, in ```/Fig1_ridge```, ```linear_simulation.R``` is used to produce the data ```linear_data_to_plot.RData```, together with ```non_linear_data_to_plot.RData``` from nonlinear data setting,  then used by ```combined_plot.R``` to produce the final figure ```combined_plot.pdf```.