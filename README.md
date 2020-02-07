# ShiftEst

## Instrumental Variable Methods using Dynamic Interventions

Repository for the Shift Estimator when using a continuous IV. This contains all the code needed
to reproduce our paper, although the data are proprietary. 

Simulations.R contains all simulation code
runOnData.R contains all analysis code -- however, these cannot be replicated without the data
doubleShiftRangeWrapper.R has the function to run the double shift estimator. To do so, specify your data, your preferred shift levels and additional arguments. This can be adapted to give other interventions through the specification of the delta vectors. 

Contact jacqueline.amauro@gmail.com with questions
