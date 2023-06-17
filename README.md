# USER GUIDE
The code for all the models developed in this project and the tests are in this repository. The dataset are put in the dataset folder and all the figures are put int the Figures folder. There is some change in figures in terms of line width, color and grids. They are adjusted for the view in the report. The use of each scripts corresponding to each test are stated below in the order presented in the report:
## Syntheic data Test
The notreal.m script is used in this part for generating synthetic data and simualting different market states. You can modify the parameters of the rand function to simulate different states.
## Parameter test
Scripts used in this test are find_epsilon.m, find_target.m, test_epsilon.m and test_target.m. The first two relates to the test of changing parameters and the last two relate to finding the best value of the parameter for the dataset.
## Computation cost test
The measure_time_horizon.m and measure_time_size.m are for the measurement of computation time for changing time horizon and portfolio size.
## Performance test
The main.m and outs.m are used for this test. The first one is for in-sample performance while the second one is for out-of-sample performance.
## Others
Those three scripts are used for the calculation of returns and cumulative returns: cal_return.m, cal_sreturn.m and cal_creturn.m. The first two take in the weight produced by the model and compute the real returns (first one for in-sample and second one for out-of-sample). The cal_creturn.m computes the cumulative returns.

If there is a need to change the dataset, apart from changing the name/address of the dataset in both the test scripts and the scripts for calculation of the returns, you also need to redefine the range of the return matrix.
