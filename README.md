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

For male analysis
- `attendance_12ponds_primarymove_males.txt`:
- `state2014_12ponds_males.txt`:
- `state2015_12ponds_males.txt`:
- `state2016_12ponds_males.txt`:
- `state2017_12ponds_males.txt`:
- `state2018_12ponds_males.txt`:
- `state2019_12ponds_males.txt`:

For female analysis
- `attendance_12ponds_primarymove_females.txt`:
- `state2014_females.txt`:
- `state2015_females.txt`:
- `state2016_females.txt`:
- `state2017_females.txt`:
- `state2018_females.txt`:
- `state2019_females.txt`:

Data matrix describing distance between 12 wetlands in study region

### Scripts 
For male analysis
- `RUN_jDist_fullmodel_M.R`:
- `FUNCTIONS_jDist_fullmodel_M.R`:
- `BOOTSTRAPS_jDist_fullmodel_M.R`:

For female analysis
- `RUN_jDist_primarymove_sec_unobs_F.R`:
- `FUNCTIONS_jDist_primarymove_sec_unobs_F.R`:
- `BOOTSTRAPS_jDist_primarymove_sec_unobs_F.R`:


