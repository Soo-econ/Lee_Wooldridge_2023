clear
set more off
*cd ""


* Lee & Wooldridge (2023, WP)
* T=6, S=4 (first treatment occur at T=4)
* Common timing intervention setup
* Estimation part only 

**************************************************************
* [1] Data Generating Process ********************************
**************************************************************

**************************************************************
* [2] Estimation**********************************************
**************************************************************

use 1.lee_wooldridge_common_data.dta, replace


xtset id year

*** Using Rolling method
	* Control group includes not-yet-treated units in this estimation.

*** Generating \dot{Y}_{it}
	bysort id: gen y_dot = y - (L1.y + L2.y + L3.y)/3 if f04
	bysort id: replace y_dot = y - (L2.y + L3.y + L4.y)/3 if f05
	bysort id: replace y_dot = y - (L3.y + L4.y + L5.y)/3 if f06

*** [ Rolling RA ] 

	teffects ra (y_dot x1 x2) (d) if f04, atet
	teffects ra (y_dot x1 x2) (d) if f05, atet
	teffects ra (y_dot x1 x2) (d) if f06, atet


*** [ Rolling IPWRA ] 
	*not-yet-treated units as well
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f04, atet
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f05, atet
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f06, atet
	
*** [ Rolling PS Matching ]		
	
	teffects psmatch (y_dot) (d x1 x2) if f04, atet
	teffects psmatch (y_dot) (d x1 x2) if f05, atet
	teffects psmatch (y_dot) (d x1 x2) if f06, atet

	* Cf.
* POLS, Wooldridge(2021):
sum x1 if d
gen x1_dm = x1 - r(mean)
sum x2 if d
gen x2_dm = x2 - r(mean)

	local second c.d#c.f04 c.d#c.f05 c.d#c.f06
	local dm c.d#c.f04#c.x1_dm c.d#c.f05#c.x1_dm c.d#c.f06#c.x1_dm c.d#c.f04#c.x2_dm c.d#c.f05#c.x2_dm c.d#c.f06#c.x2_dm 
	local dx c.d#c.x1 c.d#c.x2
	local xt_ i.year#c.x1 i.year#c.x2

reg y i.year d x1 x2 `second' `dm' `dx' `xt_', vce(cluster id)


* Callaway and Sant'Anna (2021):
	* 1) W/ never-tgreated units as control
	* First_treatment year
gen first_treat=0
	replace first_treat = 2004 if d

csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dripw) reps(0)




**************************************************************
* [3] Simulation *********************************************
**************************************************************

