clear
set more off
*cd ""

* Lee & Wooldridge (2023, WP)
* Staggered intervention setup: T=6, g={∞,4,5,6}
* Estimation part only 

**************************************************************
* [1] Data Generating Process ********************************
**************************************************************

**************************************************************
* [2] Estimation**********************************************
**************************************************************
use 2.lee_wooldridge_staggered_data.dta, replace


* Generating \dot{Y}_{igt} :
xtset id year
bysort id: gen y_44 = y - (L1.y + L2.y + L3.y)/3 if f04
bysort id: gen y_45 = y - (L2.y + L3.y + L4.y)/3 if f05
bysort id: gen y_46 = y - (L3.y + L4.y + L5.y)/3 if f06
 
bysort id: gen y_55 = y - (L1.y + L2.y + L3.y + L4.y)/4 if f05
bysort id: gen y_56 = y - (L2.y + L3.y + L4.y +L5.y)/4 if f06

bysort id: gen y_66 = y - (L1.y + L2.y + L3.y + L4.y+L5.y)/5 if f06

**Using Rolling methods
	* Control group includes not-yet-treated units in this estimation.


*** [ Rolling RA ] 

	teffects ra (y_44 x1 x2) (g4) if f04, atet
	teffects ra (y_45 x1 x2) (g4) if f05 & ~g5, atet
	teffects ra (y_46 x1 x2) (g4) if f06 & (g5 + g6 != 1), atet
	
	teffects ra (y_55 x1 x2) (g5) if f05 & ~g4, atet
	teffects ra (y_56 x1 x2) (g5) if f06 & (g4 + g6 != 1), atet
	teffects ra (y_66 x1 x2) (g6) if f06 & (g4 + g5 != 1), atet

	
*** [ Rolling IPWRA ] 

	teffects ipwra (y_44 x1 x2) (g4 x1 x2) if f04, atet
	teffects ipwra (y_45 x1 x2) (g4 x1 x2) if f05 & ~g5, atet
	teffects ipwra (y_46 x1 x2) (g4 x1 x2) if f06 & (g5 + g6 != 1), atet
	
	teffects ipwra (y_55 x1 x2) (g5 x1 x2) if f05 & ~g4, atet
	teffects ipwra (y_56 x1 x2) (g5 x1 x2) if f06 & (g4 + g6 != 1), atet
	teffects ipwra (y_66 x1 x2) (g6 x1 x2) if f06 & (g4 + g5 != 1), atet

	
*** [ Rolling PSM ]
	teffects psmatch (y_44) (g4 x1 x2) if f04, atet
	teffects psmatch (y_45) (g4 x1 x2) if f05 & ~g5, atet
	teffects psmatch (y_46) (g4 x1 x2) if f06 & (g5 + g6 != 1), atet	
	teffects psmatch (y_55) (g5 x1 x2) if f05 & ~g4, atet
	teffects psmatch (y_56) (g5 x1 x2) if f06 & (g4 + g6 != 1), atet
	teffects psmatch (y_66) (g6 x1 x2) if f06 & (g4 + g5 != 1), atet

** Cf.

* POLS, Wooldridge(2021):
	sum x1 if g4
	gen x1_dm4 = x1 - r(mean)
	sum x1 if g5
	gen x1_dm5 = x1 - r(mean)
	sum x1 if g6
	gen x1_dm6 = x1 - r(mean)
	
	sum x2 if g4
	gen x2_dm4 = x2 - r(mean)
	sum x2 if g5
	gen x2_dm5 = x2 - r(mean)
	sum x2 if g6
	gen x2_dm6 = x2 - r(mean)

	local second c.w#c.g4#c.f04 c.w#c.g4#c.f05 c.w#c.g4#c.f06
	local third c.w#c.g5#c.f05 c.w#c.g5#c.f06
	local fourth c.w#c.g6#c.f06
	local s_dm c.w#c.g4#c.f04#c.x1_dm4 c.w#c.g4#c.f05#c.x1_dm4 c.w#c.g4#c.f06#c.x1_dm4 c.w#c.g4#c.f04#c.x2_dm4 c.w#c.g4#c.f05#c.x2_dm4 c.w#c.g4#c.f06#c.x2_dm4
	local t_dm c.w#c.g5#c.f05#c.x1_dm5 c.w#c.g5#c.f06#c.x1_dm5 c.w#c.g5#c.f05#c.x2_dm5 c.w#c.g5#c.f06#c.x2_dm5
	local f_dm c.w#c.g6#c.f06#c.x1_dm6 c.w#c.g6#c.f06#c.x2_dm6
	local gx c.g4#c.x1 c.g4#c.x2 c.g5#c.x1 c.g5#c.x2 c.g6#c.x1 c.g6#c.x2
	local xt_ c.f04#c.x1 c.f05#c.x1 c.f06#c.x1 c.f04#c.x2 c.f05#c.x2 c.f06#c.x2

reg y i.year g4 g5 g6 x1 x2 `second' `third' `fourth' `s_dm' `t_dm' `f_dm' `gx' `xt_', vce(cluster id)


*** Callaway and Sant'Anna (2021):
	* 1) W/ never-tgreated units as control
	* First_treatment year
gen first_treat = 0
	replace first_treat = 2004 if g4
	replace first_treat = 2005 if g5
	replace first_treat = 2006 if g6
	
csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dripw) reps(0)



**************************************************************
* [3] Simulation *********************************************
**************************************************************
