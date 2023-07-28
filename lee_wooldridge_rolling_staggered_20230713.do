clear
set more off
log using rolling_staggered_20230713_mean_corr_ps_corr_1000, text replace

* Lee & Wooldridge (2023, WP)
* Staggered intervention setup: T=6, g={∞,4,5,6}

**************************************************************
* [1] Data Generating Process ********************************
**************************************************************

set seed 123


* # of units& replication:  n = {100, 500, 1000} , iter = 1000
global n = 1000
global iter = 1000

set obs $n

gen id = _n
expand 6
sort id

bysort id: gen year = _n + 2000
tab year, gen(f0)

* Generating two covariates (time-invariant)

*x1: asymmetrical dist.
gen x0 = rgamma(2,2)
egen x1 = mean(x0), by(id)

*x2: binary variable
gen x2 = 0.3+ rnormal(0,1) > 0
bysort id: replace x2=x2[6]
tab x2

gen c = rnormal(2,1)
bysort id: replace c = c[1]

gen u = rnormal(0,4)
gen u4 = rnormal(0,4)
gen u5 = rnormal(0,4)
gen u6 = rnormal(0,4)
sum x1 x2 c u

* Generating Treated group/cohort dummy   
* g={0,4,5,6} ; g0 -> never treated units (∞)
* g`j': first treated in 200`j' (ex. g4's initial treatment year = 2004)
* PS indicator; function of X=(x1,x2)
* 1) W/ PS correct specification
gen xgamma = -1.2+ (x1 - 4)/2 - x2
* 2) W/ PS misspecificaion
gen xgamma_m = -1.2 + (x1 - 4)/2 - x2 +((x1 - 4)^2)/2
*gen xgamma_m = -1.2+ (x1 - 4)/2 - x2 +((x1-4)^2)/5
 
* replace with "xgamma_m" for misspecification case
gen xgamma4 = xgamma
gen xgamma5 = xgamma
gen xgamma6 = xgamma


gen dnom = (1 + exp(xgamma4) + exp(xgamma5) +exp(xgamma6) )

gen p0 = 1/dnom 
gen p4 = exp(xgamma4)/dnom 
gen p5 = exp(xgamma5)/dnom 
gen p6 = exp(xgamma6)/dnom 

* Check the sum of probs. = 1
bysort id: gen unity=p0+p4+p5+p6
 
gen ru= runiform()
bysort id: replace ru=ru[1]
gen group = cond(ru <p0, 0, cond(ru<p0+p4, 4, cond(ru<p0+p4+p5, 5, 6)))
tab group,gen (g)
rename g1 g0
forval j=4(-1)2 {
local J = `j'+2
rename g`j' g`J'
}

	
* First_treatment year
gen first_treat = 0
	replace first_treat = 2004 if g4
	replace first_treat = 2005 if g5
	replace first_treat = 2006 if g6

	
** Generating potential outcomes W/ Neglected nonlinearity:
gen beta = 1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06
gen fx = (x1 - 4)/3 + x2/2 
gen hx = (x1 - 4)/2 + x2/3 
gen fx_m = (x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + (x1 -  4)*x2/4 
*gen fx_m = (x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + (x1 - 4)*x2

gen lamda_4 = 0.5*f04+0.6*f05+f06
gen lamda_5 = f05+0.5*f06
gen lamda_6 = 0.5*f06

* treatment effects : the length of exposure to the treatment.
bysort id: gen delta_t=_n
	gen te_4 = 4
	gen te_5 = 3
	gen te_6 = 2
	
forval j=4/6 {
gen hx_m_`j'=  x2*(x1-4)/3+ ((x1 - 4)^2)*`j'/3
}

forval j=4/6 {
gen te`j'= te_`j'*(0.5*(year - 200`j')+1)+(hx)*(lamda_`j') if year>=200`j'
}

*gen te`j'= te_`j'*(0.5*(year - 200`j')+1)+(hx_m_`j')*(lamda_`j') if year>=200`j'




** Potential outcomes with common trends
*1) Control State

gen y0 = delta_t + (fx)*(beta) +c +u
*gen y0 = delta_t + (fx_m)*(beta) +c +u

*2) Treated State
forval j=4/6 {
gen y`j'= y0
replace y`j'= y0 +te`j'-u+u`j' if year>= 200`j'
}

** Generating by hand
*gen y0 = delta_t+ ((x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + (x1 - 4)*x2)*(1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06)+ c + u
*gen y4 = y0
*replace y4 = y0 +te*(0.5*(year - first_treat) +1) +(x2*(x1-4)/3+ ((x1 - 4)^2)*4/3)*(0.5*f04+0.6*f05+f06)-u+u4 if year >= 2004
*gen y5 = y0
*replace y5 = y0 +te*(0.5*( year - first_treat) +1)  +(x2*(x1-4)/3+ ((x1 - 4)^2)*5/3)*(f05+0.5*f06)-u+u5 if year >= 2005
*gen y6 = y0
*replace y6 = y0 +te+(((x1 - 4)^2)*2 + x2*(x1-4)/3)*0.5*f06 -u+u6 if year == 2006



* Observed outcome:
gen y = g0*y0 + g4*y4 + g5*y5 + g6*y6


* Generate time-varying treatment indicator for staggered intervention:
gen w = g4*(f04 + f05 + f06) + g5*(f05 + f06) + g6*f06


*Sample att
gen att_44 = y4-y0 if year==2004 & g4
gen att_45 = y4-y0 if year==2005 & g4 & ~g5
gen att_46 = y4-y0 if year==2006 & g4 & (g5 + g6 != 1)
gen att_55 = y5-y0 if year==2005 & g5 & ~g4
gen att_56 = y5-y0 if year==2006 & g5 & (g4 + g6 != 1)
gen att_66 = y6-y0 if year==2006 & g6 & (g4 + g5 != 1) 

			sum att_4* if g4
			sum att_5* if g5
			sum att_6* if g6




**************************************************************
* [2] Estimation**********************************************
**************************************************************


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


* POLS:
xtset id year
reg y i.year g4 g5 g6 x1 x2 `second' `third' `fourth' `s_dm' `t_dm' `f_dm' `gx' `xt_', vce(cluster id)


*** Callaway and Sant'Anna:
	* 1) W/ only never-tgreated units as control
csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dripw) reps(0)


* IPWRA using all control periods:
bysort id: gen y_44 = y - (L1.y + L2.y + L3.y)/3 if f04
bysort id: gen y_45 = y - (L2.y + L3.y + L4.y)/3 if f05
bysort id: gen y_46 = y - (L3.y + L4.y + L5.y)/3 if f06
 
bysort id: gen y_55 = y - (L1.y + L2.y + L3.y + L4.y)/4 if f05
bysort id: gen y_56 = y - (L2.y + L3.y + L4.y +L5.y)/4 if f06

bysort id: gen y_66 = y - (L1.y + L2.y + L3.y + L4.y+L5.y)/5 if f06

**Use an alternative rolling method

***[ Rolling RA ] 
	* W/ not-yet-treated units as well
	teffects ra (y_44 x1 x2) (g4) if f04, atet
	teffects ra (y_45 x1 x2) (g4) if f05 & ~g5, atet
	teffects ra (y_46 x1 x2) (g4) if f06 & (g5 + g6 != 1), atet
	
	teffects ra (y_55 x1 x2) (g5) if f05 & ~g4, atet
	teffects ra (y_56 x1 x2) (g5) if f06 & (g4 + g6 != 1), atet
	teffects ra (y_66 x1 x2) (g6) if f06 & (g4 + g5 != 1), atet

	
*** [ Rolling IPWRA ] 

	*not-yet-treated units as well
	teffects ipwra (y_44 x1 x2) (g4 x1 x2) if f04, atet
	teffects ipwra (y_45 x1 x2) (g4 x1 x2) if f05 & ~g5, atet
	teffects ipwra (y_46 x1 x2) (g4 x1 x2) if f06 & (g5 + g6 != 1), atet
	
	teffects ipwra (y_55 x1 x2) (g5 x1 x2) if f05 & ~g4, atet
	teffects ipwra (y_56 x1 x2) (g5 x1 x2) if f06 & (g4 + g6 != 1), atet
	teffects ipwra (y_66 x1 x2) (g6 x1 x2) if f06 & (g4 + g5 != 1), atet

	
*** [ Rolling IPWRA ]
	
	teffects psmatch (y_44) (g4 x1 x2) if f04, atet
	teffects psmatch (y_45) (g4 x1 x2) if f05 & ~g5, atet
	teffects psmatch (y_46) (g4 x1 x2) if f06 & (g5 + g6 != 1), atet	
	teffects psmatch (y_55) (g5 x1 x2) if f05 & ~g4, atet
	teffects psmatch (y_56) (g5 x1 x2) if f06 & (g4 + g6 != 1), atet
	teffects psmatch (y_66) (g6 x1 x2) if f06 & (g4 + g5 != 1), atet


**************************************************************
* [3] Simulation *********************************************
**************************************************************

capture program drop did_T6

program did_T6, rclass
drop _all

set obs $n

gen id = _n
expand 6
sort id

bysort id: gen year = _n + 2000
tab year, gen(f0)

* Generating two covariates (time-invariant)

*x1: asymmetrical dist.
gen x0 = rgamma(2,2)
egen x1 = mean(x0), by(id)

*x2: binary variable
gen x2 = 0.3+ rnormal(0,1) > 0
bysort id: replace x2=x2[6]
tab x2

gen c = rnormal(2,1)
bysort id: replace c = c[1]

gen u = rnormal(0,4)
gen u4 = rnormal(0,4)
gen u5 = rnormal(0,4)
gen u6 = rnormal(0,4)
sum x1 x2 c u




* Generating Treated group/cohort dummy   
* g={0,4,5,6} ; g0 -> never treated units (∞)
* g`j': first treated in 200`j' (ex. g4's initial treatment year = 2004)
* PS indicator; function of X=(x1,x2)
* 1) W/ PS correct specification
gen xgamma = -1.2+ (x1 - 4)/2 - x2
* 2) W/ PS misspecificaion
gen xgamma_m = -1.2 + (x1 - 4)/2 - x2 +((x1 - 4)^2)/2


* replace with xgamma_m for misspecification case
gen xgamma4 = xgamma
gen xgamma5 = xgamma
gen xgamma6 = xgamma

gen dnom = (1 + exp(xgamma4) + exp(xgamma5) +exp(xgamma6) )

gen p0 = 1/dnom 
gen p4 = exp(xgamma4)/dnom 
gen p5 = exp(xgamma5)/dnom 
gen p6 = exp(xgamma6)/dnom 

* Check the sum of probs. = 1
bysort id: gen unity=p0+p4+p5+p6
 
gen ru= runiform()
bysort id: replace ru=ru[1]
gen group = cond(ru <p0, 0, cond(ru<p0+p4, 4, cond(ru<p0+p4+p5, 5, 6)))
tab group,gen (g)
rename g1 g0
forval j=4(-1)2 {
local J = `j'+2
rename g`j' g`J'
}




	*************************************
		sum g0
			return scalar g0_p = r(mean)
			sum g4
			return scalar g4_p = r(mean)
			sum g5
			return scalar g5_p = r(mean)
			sum g6
			return scalar g6_p = r(mean)
	*************************************


	
* First_treatment year
gen first_treat = 0
	replace first_treat = 2004 if g4
	replace first_treat = 2005 if g5
	replace first_treat = 2006 if g6
	


** Generate potential outcomes
** 1) W/ Neglected nonlinearity (misspecificaion: fx_m/ hx_m)
gen beta = 1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06
gen fx = (x1 - 4)/3 + x2/2 
gen hx = (x1 - 4)/2 + x2/3 
gen fx_m = (x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + (x1 -  4)*x2/4 
*gen fx_m = (x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + (x1 - 4)*x2

gen lamda_4 = 0.5*f04+0.6*f05+f06
gen lamda_5 = f05+0.5*f06
gen lamda_6 = 0.5*f06

* treatment effect (att): is the length of exposure to the treatment.
bysort id: gen delta_t=_n
	gen te_4 = 4
	gen te_5 = 3
	gen te_6 = 2
	
forval j=4/6 {
gen hx_m_`j'=  x2*(x1-4)/3+ ((x1 - 4)^2)*`j'/3
}

forval j=4/6 {
gen te`j'= te_`j'*(0.5*(year - 200`j')+1)+(hx)*(lamda_`j') if year>=200`j'
}

*gen te`j'= te_`j'*(0.5*(year - 200`j')+1)+(hx_m_`j')*(lamda_`j') if year>=200`j'


*** Potential outcomes
		* 1) control state
		gen y0 = delta_t + (fx)*(beta) +c +u
		*gen y0 = delta_t + (fx_m)*(beta) +c +u
		* 2) treated state
		forval j=4/6 {
		gen y`j'= y0
		replace y`j'= y0 +te`j'-u+u`j' if year>= 200`j'
		}

*** Observed outcome:
gen y = g0*y0 + g4*y4 + g5*y5 + g6*y6


*** Generate time-varying treatment indicator for staggered intervention:
gen w = g4*(f04 + f05 + f06) + g5*(f05 + f06) + g6*f06


*** Sample ATT 
gen att_44 = y4-y0 if year==2004 & g4
gen att_45 = y4-y0 if year==2005 & g4 & ~g5
gen att_46 = y4-y0 if year==2006 & g4 & (g5 + g6 != 1)
gen att_55 = y5-y0 if year==2005 & g5 & ~g4
gen att_56 = y5-y0 if year==2006 & g5 & (g4 + g6 != 1)
gen att_66 = y6-y0 if year==2006 & g6 & (g4 + g5 != 1) 


********************************************************
			sum att_44 if g4
			return scalar att_44 = r(mean)
			sum att_45 if g4
			return scalar att_45 = r(mean)
			sum att_46 if g4
			return scalar att_46 = r(mean)
	
	
			sum att_55 if g5
			return scalar att_55 = r(mean)
			sum att_56 if g5
			return scalar att_56 = r(mean)
	
			sum att_66 if g6
			return scalar att_66 = r(mean)
			
********************************************************

		
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

local second c.g4#c.f04 c.g4#c.f05 c.g4#c.f06
local third c.g5#c.f05 c.g5#c.f06
local fourth c.g6#c.f06
local s_dm c.g4#c.f04#c.x1_dm4 c.g4#c.f05#c.x1_dm4 c.g4#c.f06#c.x1_dm4 c.g4#c.f04#c.x2_dm4 c.g4#c.f05#c.x2_dm4 c.g4#c.f06#c.x2_dm4
local t_dm c.g5#c.f05#c.x1_dm5 c.g5#c.f06#c.x1_dm5 c.g5#c.f05#c.x2_dm5 c.g5#c.f06#c.x2_dm5
local f_dm c.g6#c.f06#c.x1_dm6 c.g6#c.f06#c.x2_dm6
local gx c.g4#c.x1 c.g4#c.x2 c.g5#c.x1 c.g5#c.x2 c.g6#c.x1 c.g6#c.x2
local xt_ c.f04#c.x1 c.f05#c.x1 c.f06#c.x1 c.f04#c.x2 c.f05#c.x2 c.f06#c.x2


*plot d x1_dm

* POLS:
xtset id year
reg y i.year g4 g5 g6 x1 x2 `second' `third' `fourth' `s_dm' `t_dm' `f_dm' `gx' `xt_', vce(cluster id)

	return scalar att_44_pols = _b[c.g4#c.f04]
	return scalar att_45_pols = _b[c.g4#c.f05]
	return scalar att_46_pols = _b[c.g4#c.f06]
	return scalar att_55_pols = _b[c.g5#c.f05]
	return scalar att_56_pols = _b[c.g5#c.f06]
	return scalar att_66_pols = _b[c.g6#c.f06]
	return scalar rsq = e(r2)
	

		
* Callaway and Sant'Anna:

	csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dripw) reps(0)

	return scalar att_44_cs = _b[g2004:t_2003_2004]
	return scalar att_45_cs = _b[g2004:t_2003_2005]
	return scalar att_46_cs = _b[g2004:t_2003_2006]
	return scalar att_55_cs = _b[g2005:t_2004_2005]
	return scalar att_56_cs = _b[g2005:t_2004_2006]
	return scalar att_66_cs = _b[g2006:t_2005_2006]


* Rolling Method/ Generating Y_dot 

gen y_44 = y - (L1.y + L2.y + L3.y)/3 if f04
gen y_45 = y - (L2.y + L3.y + L4.y)/3 if f05
gen y_46 = y - (L3.y + L4.y + L5.y)/3 if f06

gen y_55 = y - (L1.y + L2.y + L3.y + L4.y)/4 if f05
gen y_56 = y - (L2.y + L3.y + L4.y +L5.y)/4 if f06

gen y_66 = y - (L1.y + L2.y + L3.y + L4.y+L5.y)/5 if f06


*** [ Rolling RA ]

	* W/ not-yet-treated as controls

	teffects ra (y_44 x1 x2) (g4) if f04, atet
	return scalar att_44_ra = _b[ATET:r1vs0.g4]
	teffects ra (y_45 x1 x2) (g4) if f05 & ~g5, atet
	return scalar att_45_ra = _b[ATET:r1vs0.g4]
	teffects ra (y_46 x1 x2) (g4) if f06 & (g5 + g6 != 1), atet
	return scalar att_46_ra = _b[ATET:r1vs0.g4]
	
	teffects ra (y_55 x1 x2) (g5) if f05 & ~g4, atet
	return scalar att_55_ra = _b[ATET:r1vs0.g5]
	teffects ra (y_56 x1 x2) (g5) if f06 & (g4 + g6 != 1), atet
	return scalar att_56_ra = _b[ATET:r1vs0.g5]
	
	teffects ra (y_66 x1 x2) (g6) if f06 & (g4 + g5 != 1), atet
	return scalar att_66_ra = _b[ATET:r1vs0.g6]	
				

*** [ Rolling IPWRA ]
		
	* W/ not-yet-treated as controls
	teffects ipwra (y_44 x1 x2) (g4 x1 x2) if f04, atet
	return scalar att_44_ipwra = _b[ATET:r1vs0.g4]
	teffects ipwra (y_45 x1 x2) (g4 x1 x2) if f05 & ~g5, atet
	return scalar att_45_ipwra = _b[ATET:r1vs0.g4]
	teffects ipwra (y_46 x1 x2) (g4 x1 x2) if f06 & (g5 + g6 != 1), atet
	return scalar att_46_ipwra = _b[ATET:r1vs0.g4]
	
	teffects ipwra (y_55 x1 x2) (g5 x1 x2) if f05 & ~g4, atet
	return scalar att_55_ipwra = _b[ATET:r1vs0.g5]
	teffects ipwra (y_56 x1 x2) (g5 x1 x2) if f06 & (g4 + g6 != 1), atet
	return scalar att_56_ipwra = _b[ATET:r1vs0.g5]
	
	teffects ipwra (y_66 x1 x2) (g6 x1 x2) if f06 & (g4 + g5 != 1), atet
	return scalar att_66_ipwra = _b[ATET:r1vs0.g6]	
		
	**** W/ never-treated (using multinomial logit pr.)
*	teffects ipwra (y_44 x1 x2) (group x1 x2) if f04, atet
*	return scalar att_44_ipwrac = _b[ATET:r4vs0.group]
*	teffects ipwra (y_45 x1 x2) (group x1 x2) if f05, atet
*	return scalar att_45_ipwrac = _b[ATET:r4vs0.group]
*	teffects ipwra (y_46 x1 x2) (group x1 x2) if f06, atet
*	return scalar att_46_ipwrac = _b[ATET:r4vs0.group]
*	teffects ipwra (y_55 x1 x2) (group x1 x2) if f05, atet
*	return scalar att_55_ipwrac = _b[ATET:r5vs0.group]
*	teffects ipwra (y_56 x1 x2) (group x1 x2) if f06, atet
*	return scalar att_56_ipwrac = _b[ATET:r5vs0.group]
	
*	teffects ipwra (y_66 x1 x2) (group x1 x2) if f06, atet
*	return scalar att_66_ipwrac = _b[ATET:r6vs0.group]	

*** [ Rolling PS Matching ]		

	teffects psmatch (y_44) (g4 x1 x2) if f04, atet
	return scalar att_44_psmatch = _b[ATET:r1vs0.g4]
	teffects psmatch (y_45) (g4 x1 x2) if f05 & ~g5, atet
	return scalar att_45_psmatch = _b[ATET:r1vs0.g4]
	teffects psmatch (y_46) (g4 x1 x2) if f06 & (g5 + g6 != 1), atet
	return scalar att_46_psmatch = _b[ATET:r1vs0.g4]
	
	teffects psmatch (y_55) (g5 x1 x2) if f05 & ~g4, atet
	return scalar att_55_psmatch = _b[ATET:r1vs0.g5]
	teffects psmatch (y_56) (g5 x1 x2) if f06 & (g4 + g6 != 1), atet
	return scalar att_56_psmatch = _b[ATET:r1vs0.g5]
	teffects psmatch (y_66) (g6 x1 x2) if f06 & (g4 + g5 != 1), atet
	return scalar att_66_psmatch = _b[ATET:r1vs0.g6]




end

set seed 123



simulate r(att_44) r(att_44_pols) r(att_44_ra) r(att_44_cs) r(att_44_ipwra) r(att_44_psmatch)  ///
		 r(att_45) r(att_45_pols) r(att_45_ra) r(att_45_cs) r(att_45_ipwra) r(att_45_psmatch)  ///
		 r(att_46) r(att_46_pols) r(att_46_ra) r(att_46_cs) r(att_46_ipwra) r(att_46_psmatch)  ///
		 r(att_55) r(att_55_pols) r(att_55_ra) r(att_55_cs) r(att_55_ipwra) r(att_55_psmatch) ///
		 r(att_56) r(att_56_pols) r(att_56_ra) r(att_56_cs) r(att_56_ipwra) r(att_56_psmatch)  ///
		 r(att_66) r(att_66_pols) r(att_66_ra) r(att_66_cs) r(att_66_ipwra) r(att_66_psmatch)  ///
		 r(rsq) r(g0_p) r(g4_p) r(g5_p) r(g6_p),reps($iter): did_T6
sum, sep(6)


log close
