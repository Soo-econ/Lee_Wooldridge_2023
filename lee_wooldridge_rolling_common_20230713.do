clear
set more off
log using rolling_common_20230713_mean_miss_ps_corr_500, text replace


* Lee & Wooldridge (2023, WP)
* T=6, S=4 (first treatment occur at T=4)
* Common timing intervention setup

**************************************************************
* [1] Data Generating Process ********************************
**************************************************************

set seed 123

* # of units& replication:  n = {100, 500, 1000} , iter = 1000
global n = 500
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

*x2: binary variable (bernoulli (0.6))
gen x2 = rbinomial(1,0.6)
*gen x2 = 0.3+ rnormal(0,1) > 0
bysort id: replace x2 = x2[6]
sum x2

gen c = rnormal(2,1)
bysort id: replace c = c[1]

gen u = rnormal(0,4)
gen u1 = rnormal(0,4)
gen u4 = rnormal(0,4)
gen u5 = rnormal(0,4)
gen u6 = rnormal(0,4)
sum x1 x2 c u

* Generating Treatment indicator

* 1) W/ PS correct specification
gen z = -1.2+ (x1 - 4)/2 - x2
* 2) W/ PS misspecificaion
*gen zm = -1.2 + (x1 - 4)/2 - x2 +((x1 - 4)^2)/2

gen pr = 1/(1 + exp(-z))
sum pr
gen d= rbinomial(1,pr)
bysort id: replace d= d[1] 
tab d 
	

	
** Generating potential outcomes W/ Neglected nonlinearity:
gen beta = 1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06

**correctly specified: fx,hx
gen fx = (x1 - 4)/3 + x2/2
gen hx = (x1 - 4)/2 + x2/3 
** W/ neglected nonlinearity: fx_m, gx_m

gen fx_m = (x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + (x1 -  4)*x2/4 
gen hx_m = (x1 - 4)/2 + x2/3 +((x1 - 4)^2)/4 + (x1 - 4)*x2/2
gen lamda = 0.5*f04+0.6*f05+f06

bysort id: gen delta_t=_n

* treatment effects : the length of exposure to the treatment.
gen te_4 = 3
gen te_5 = 4.5
gen te_6 = 5.5
gen te =0

forval j= 4/6 {
replace te = te_`j' +(hx_m)*(lamda) if year==200`j'
}


** Potential outcomes with common trends
*1) Control State
gen y0 = delta_t + (fx_m)*(beta) +c +u


*2) Treated State
gen y1= y0
replace y1= y0 +te -u+u1 if year>= 2004

** Generating by hand
*gen y0 = delta_t+ ((x1 - 4)/3 + x2/2 )*(1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06)+c + u
*gen y1 = y0
*replace y1 = y0 +att+((x1 - 4)/2+x2/3)*(0.5*f04+0.6*f05+f06)-u+u2 if year >= 2004


* Observed outcome:
gen y = (1 - d)*y0 + d*y1


*Sample att
gen te_i = y1-y0
sum te_i if d & f04
sum te_i if d & f05
sum te_i if d & f06
				

				
				

**************************************************************
* [2] Estimation**********************************************
**************************************************************

sum x1 if d
gen x1_dm = x1 - r(mean)
sum x2 if d
gen x2_dm = x2 - r(mean)

	local second c.d#c.f04 c.d#c.f05 c.d#c.f06
	local dm c.d#c.f04#c.x1_dm c.d#c.f05#c.x1_dm c.d#c.f06#c.x1_dm c.d#c.f04#c.x2_dm c.d#c.f05#c.x2_dm c.d#c.f06#c.x2_dm 
	local dx c.d#c.x1 c.d#c.x2
	local xt_ i.year#c.x1 i.year#c.x2

* POLS:
xtset id year
reg y i.year d x1 x2 `second' `dm' `dx' `xt_', vce(cluster id)


*** Callaway and Sant'Anna:
	* 1) W/ only never-tgreated units as control
	* First_treatment year
	gen first_treat=0
	replace first_treat = 2004 if d

	csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dripw) reps(0)


*** Using Rolling method
*** Generating y_dot
	bysort id: gen y_dot = y - (L1.y + L2.y + L3.y)/3 if f04
	bysort id: replace y_dot = y - (L2.y + L3.y + L4.y)/3 if f05
	bysort id: replace y_dot = y - (L3.y + L4.y + L5.y)/3 if f06


*** [ Rolling RA ] 
	* W/ not-yet-treated units as well

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

*x2: binary variable (bernoulli (0.6))
gen x2 = rbinomial(1,0.6)
*gen x2 = 0.3+ rnormal(0,1) > 0
bysort id: replace x2 = x2[6]
sum x2

gen c = rnormal(2,1)
bysort id: replace c = c[1]

gen u = rnormal(0,4)
gen u1 = rnormal(0,4)
gen u4 = rnormal(0,4)
gen u5 = rnormal(0,4)
gen u6 = rnormal(0,4)
sum x1 x2 c u

* Generating Treatment indicator

* 1) W/ PS correct specification
gen z = -1.2+ (x1 - 4)/2 - x2
* 2) W/ PS misspecificaion
*gen zm = -1.2 + (x1 - 4)/2 - x2 +((x1 - 4)^2)/2

gen pr = 1/(1 + exp(-z))
sum pr
gen d= rbinomial(1,pr)
bysort id: replace d= d[1] 
tab d 
	
sum d if f04
return scalar d_p = r(mean)


	
** Generating potential outcomes W/ Neglected nonlinearity:

gen beta = 1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06
**correctly specified: fx,hx
gen fx = (x1 - 4)/3 + x2/2 
gen hx = (x1 - 4)/2 + x2/3

** W/ neglected nonlinearity: fx_m, gx_m
gen fx_m = (x1 - 4)/3 + x2/2 +((x1 - 4)^2)/3 + x2*(x1 -  4)/4 
gen hx_m = (x1 - 4)/5 + x2/3 +((x1 - 4)^2)/4 + (x1 - 4)*x2/2

gen lamda = 0.5*f04+0.6*f05+f06

bysort id: gen delta_t=_n

* treatment effects : the length of exposure to the treatment.
gen te_4 = 3
gen te_5 = 4.5
gen te_6 = 5.5
gen te =0

forval j= 4/6 {
replace te = te_`j' +(hx_m)*(lamda) if year==200`j'
}


** Potential outcomes with common trends
*1) Control State
*gen y0 = delta_t + (fx)*(beta) +c +u
gen y0 = delta_t + (fx_m)*(beta) + c +u

*2) Treated State
gen y1= y0
replace y1= y0 +te -u+u1 if year>= 2004

** Generating by hand
*gen y0 = delta_t+ ((x1 - 4)/2 + x2/3 )*(1*f01+1.5*f02+0.8*f03+1.5*f04+2*f05+2.5*f06)+c + u
*gen y1 = y0
*replace y1 = y0 +att+((x1 - 4)/3+x2/2)*(0.5*f04+0.6*f05+f06)-u+u2 if year >= 2004


* Observed outcome:
gen y = (1 - d)*y0 + d*y1


*Sample att
gen te_i = y1-y0

********************************************************
	sum te_i if d & f04
	return scalar att_4 = r(mean)
	
	sum te_i if d & f05
	return scalar att_5 = r(mean)

	sum te_i if d & f06
	return scalar att_6 = r(mean)

			
********************************************************


sum x1 if d
gen x1_dm = x1 - r(mean)
sum x2 if d
gen x2_dm = x2 - r(mean)

	local second c.d#c.f04 c.d#c.f05 c.d#c.f06
	local dm c.d#c.f04#c.x1_dm c.d#c.f05#c.x1_dm c.d#c.f06#c.x1_dm c.d#c.f04#c.x2_dm c.d#c.f05#c.x2_dm c.d#c.f06#c.x2_dm
	local dx c.d#c.x1 c.d#c.x2
	local xt_ i.year#c.x1 i.year#c.x2	
	
*** POLS:
xtset id year
reg y i.year d x1 x2 `second' `dm' `dx' `xt_', vce(cluster id)

	return scalar att_4_pols = _b[c.d#c.f04]
	return scalar att_5_pols = _b[c.d#c.f05]
	return scalar att_6_pols = _b[c.d#c.f06]
	return scalar rsq = e(r2)

*** Callaway and Sant'Anna:
	* 1) W/ only never-tgreated units as control
	* First_treatment year
	gen first_treat=0
	replace first_treat = 2004 if d

	csdid y x1 x2, ivar(id) time(year) gvar(first_treat) method(dripw) reps(0)

	return scalar att_4_cs = _b[g2004:t_2003_2004]
	return scalar att_5_cs = _b[g2004:t_2003_2005]
	return scalar att_6_cs = _b[g2004:t_2003_2006]	
	
	
*** Rolling method 
*** Generating ydot
	bysort id: gen y_dot = y - (L1.y + L2.y + L3.y)/3 if f04
	bysort id: replace y_dot = y - (L2.y + L3.y + L4.y)/3 if f05
	bysort id: replace y_dot = y - (L3.y + L4.y + L5.y)/3 if f06
	
*** [ Rolling RA with correct functional form ]

	teffects ra (y_dot x1 x2 c.x1#c.x1 c.x1#c.x2) (d) if f04, atet
	return scalar att_4_ra_full = _b[ATET:r1vs0.d]
	teffects ra (y_dot x1 x2 c.x1#c.x1 c.x1#c.x2) (d) if f05, atet
	return scalar att_5_ra_full = _b[ATET:r1vs0.d]
	teffects ra (y_dot x1 x2 c.x1#c.x1 c.x1#c.x2) (d) if f06, atet
	return scalar att_6_ra_full = _b[ATET:r1vs0.d]


*** [ Rolling IPWRA ]
		

	teffects ipwra (y_dot x1 x2) (d x1 x2) if f04, atet
	return scalar att_4_ipwra = _b[ATET:r1vs0.d]
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f05, atet
	return scalar att_5_ipwra = _b[ATET:r1vs0.d]
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f06, atet
	return scalar att_6_ipwra = _b[ATET:r1vs0.d]
		
*** [ Rolling PS Matching ]		

	teffects psmatch (y_dot) (d x1 x2) if f04, atet
	return scalar att_4_psmatch = _b[ATET:r1vs0.d]
	teffects psmatch (y_dot) (d x1 x2) if f05, atet
	return scalar att_5_psmatch = _b[ATET:r1vs0.d]
	teffects psmatch (y_dot) (d x1 x2) if f06, atet
	return scalar att_6_psmatch = _b[ATET:r1vs0.d]


end

set seed 123



simulate r(att_4) r(att_4_pols) r(att_4_ra_full) r(att_4_cs) r(att_4_ipwra) r(att_4_psmatch) ///
	r(att_5) r(att_5_pols) r(att_5_ra_full) r(att_5_cs) r(att_5_ipwra) r(att_5_psmatch) ///
	r(att_6) r(att_6_pols) r(att_6_ra_full) r(att_6_cs) r(att_6_ipwra) r(att_6_psmatch) ///
	r(rsq) r(d_p), reps($iter): did_T6
sum, sep(6)

log close
