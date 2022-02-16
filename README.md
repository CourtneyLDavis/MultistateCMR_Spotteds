# Multistate mark-recapture model to estimate sex-specific dispersal rates and distances for a wetland breeding amphibian metapopulation
In Review 

Courtney L. Davis*, David J. Muñoz, Staci M. Amburgey, Carli R. Dinsmore, Eric W. Teitsworth, David A. W. Miller

*Please contact the first author with any questions about the code or publication: Courtney L. Davis (cld74@cornell.edu)


## Abstract
1. How animals move across space and time in heterogeneous landscapes has important implications for the conservation of imperiled and dispersal-limited taxa, such as amphibians. Wetland-breeding amphibians are thought to exhibit strong site fidelity, but there is growing evidence that many species are more vagile than previously assumed. Few studies have quantified breeding dispersal probabilities or distances while also accounting for observational uncertainty, which means that resultant estimates could be biased or confounded with sampling intensity and other demographic rates (e.g., site and time specific detection and survival probabilities).

2. Here, we analyze data from a 6-year capture-mark-recapture study on adult spotted salamanders (Ambystoma maculatum) conducted at 12 wetlands in central Pennsylvania to estimate population dynamics and metapopulation structure at multiple scales. 

3. We used a multistate, hidden Markov estimator to quantify sex-specific site fidelity and breeding dispersal as a function of Euclidean distance between wetlands, while accounting for imperfect detection. We estimated short-time scale movements (i.e., those that occur within a single breeding season) and longer-time scale movements (i.e., those that occur among breeding seasons) to determine if dispersal rates and distances differed. 

4. We found that inter-annual site fidelity of males varied among wetlands and was positively associated with population density. Females exhibited higher inter-annual site fidelity and dispersed further than males between breeding seasons. Within breeding seasons, we found that up to 6% of males dispersed to a new wetland each day. 

## Repository Directory

### Data
This contains all of the data files necessary to run the spotted salamander analysis for males and females. It includes the following:

For male analysis where movement was modeled between and within years as a function of distance between sites
- `attendance_12ponds_primarymove_males.txt`: data file describing the pond (sites 1-12) where each individual was first captured in each year of the study.
- `state2014_12ponds_males.txt`: data file describing the pond (sites 7-9) where each individual was captured on each trap day (n = 6 days) in 2014. Note that only 3 ponds were trapped in 2014.
- `state2015_12ponds_males.txt`: data file describing the pond (sites 1-12) where each individual was captured on each trap day (n = 15 days) in 2015
- `state2016_12ponds_males.txt`: data file describing the pond (sites 1-12) where each individual was captured on each trap day (n = 17 days) in 2016
- `state2017_12ponds_males.txt`: data file describing the pond (sites 1-12) where each individual was captured on each trap day (n = 15 days) of 2017
- `state2018_12ponds_males.txt`: data file describing the pond (sites 1-12) where each individual was captured on each trap day (n = 21 days) of 2018
- `state2019_12ponds_males.txt`: data file describing the pond (sites 1-12) where each individual was captured on each trap day (n = 10 days) of 2019

For female analysis where movement was only modeled between years as a function of distance between sites
- `attendance_12ponds_primarymove_females.txt`: data file descibing the pond (sites 1-12) where each individual was first captured in each year of the study.
- `state2014_females.txt`: data file describing whether an individual was captured on each trap day (i.e., 1 if captured, 0 if not) in 2014
- `state2015_females.txt`: data file describing whether an individual was captured on each trap day (i.e., 1 if captured, 0 if not) in 2015
- `state2016_females.txt`: data file describing whether an individual was captured on each trap day (i.e., 1 if captured, 0 if not) in 2016
- `state2017_females.txt`: data file describing whether an individual was captured on each trap day (i.e., 1 if captured, 0 if not) in 2017
- `state2018_females.txt`: data file describing whether an individual was captured on each trap day (i.e., 1 if captured, 0 if not) in 2018
- `state2019_females.txt`: data file describing whether an individual was captured on each trap day (i.e., 1 if captured, 0 if not) in 2019

Distance between sites
- `DistanceMatrix_AMMACMR.csv`: data file describing the distance between the 12 ponds in our study region


### Scripts 
For male analysis
- `RUN_jDist_fullmodel_M.R`: R script that defines the data, sets initial values, calls `FUNCTIONS_jDist_fullmodel_M.R` to fit the model, and calls `BOOTSTRAPS_jDist_fullmodel_M.R` to bootstrap results.
- `FUNCTIONS_jDist_fullmodel_M.R`: R script that unpacks the parameter vector and defines the model likelihood. 
- `BOOTSTRAPS_jDist_fullmodel_M.R`: R script that performs the bootstrap to calculate bootstrap standard errors for all defined parameters. 

For female analysis
- `RUN_jDist_primarymove_sec_unobs_F.R`: R script that defines the data, sets initial values, calls `FUNCTIONS_jDist_primarymove_sec_unobs_F.R` to fit the model, and calls `BOOTSTRAPS_jDist_primarymove_sec_unobs_F.R` to bootstrap results.
- `FUNCTIONS_jDist_primarymove_sec_unobs_F.R`: R script that unpacks the parameter vector and defines the model likelihood. 
- `BOOTSTRAPS_jDist_primarymove_sec_unobs_F.R`: R script that performs the bootstrap to calculate bootstrap standard errors for all defined parameters. 


