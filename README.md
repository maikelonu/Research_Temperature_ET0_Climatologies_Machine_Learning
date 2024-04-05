# CRCLIMtemp-TEC

CRCLIMtemp is an program intended for the generation of temperature/evapotranspiration climatologies using several interpolation methods

## Installation

Use Git to clone and install the program

## Usage

Run PAPER_MMENDEZ_2018_10.R script accordingly

## Contributing

Maikel Mendez Morales. Escuela de Ingeniería en Construcción, Instituto Tecnológico de Costa Rica. email: mamendez@itcr.ac.cr

Luis Alexander Calvo Valverde. Escuela de Ingeniería en Computación, Instituto Tecnológico de Costa Rica. email: lcalvo@itcr.ac.cr

## Publications

Comparison performance of machine learning and geostatistical methods for the interpolation of monthly air temperature over Costa Rica

IOP Conference Series: Earth and Environmental Science

https://iopscience.iop.org/article/10.1088/1755-1315/432/1/012011/meta

https://doi.org/10.1088/1755-1315/432/1/012011

![alt test](/mansc_721.png)

Graphical Abstract

![alt test](/QGIS_Temp_ET0.png)

Abstract: 

The performance of three machine learning (ML) methods; cubist regression (CR), random forest (RF) and generalized additive model using splines (GAM) in generating monthly air temperature grids over Costa Rica was evaluated against two heavily used geostatistical methods; ordinary kriging (OK) and kriging with external drift (KED). The skill of the interpolation methods was evaluated using a 10-fold cross-validation technique; selecting the root-mean square error (RMSE), the mean absolute error (MAE) and the Pearson correlation-coefficient (R) as agreement metrics. To this purpose, data from an irregularlydistributed observational-network comprised by 73 weather-stations were selected for the period 1950-1987. Several spatial fields derived from a high-resolution digital elevation model (DEM) were tested as covariants. Results from the 10-fold cross-validation test show that CR
yielded the best individual performance followed by KED, whereas GAM performed worst. Elevation on the other hand, was the only covariant ultimately incorporated in the interpolation process, since the remaining spatial fields exhibited poor correlation with temperature or resulted in data redundancy. While the quantitative and qualitative evaluation of CR and KED
can be said to be comparable, CR is considered the best approach since the method is unaffected by assumptions on data normality and homoscedasticity.

Results

![alt test](/compiled_temperature_1951-1990.png)

![alt test](/compiled_et0_1951-1990.png)

![alt test](/FR18_01.png)

## License

[MIT](https://choosealicense.com/licenses/mit/)

