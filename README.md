## Lee and Wooldridge 2023
_"A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data"_  
	Available on SSRN https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518

## How to Apply Rolling Methods to your Panal Data
## 1. Common Timing Case

__Basic set up__

1) Time Periods: $t \in {1,...,T}$

2) Time dummies: $f2_t, ...,fT_t$ , i.e., $f4 =1$ if $t=4$

3) The first Intervention occurs at $S$, $1 < S \leq T$

4) There are two observable covariates $X=(X_1, X_2)$

### Procedure 3.1
__Step.1__ For a given time period $t = S, \ldots, T$ and each unit $i$, compute
```math
    \dot{Y}_{it} \equiv Y_{it}-\frac{1}{S-1} \sum_{q=1}^{S-1} Y_{iq}
```

__Step.2__ Using all of units, apply standard TE methods - such as linear RA, IPW, IPWRA, PS matching - to the cross section
```math
\{ ( \dot{Y}_{it}, D_i, \mathbf{X}_i) \ , \ i \ = \ 1, \ldots, N ;  t= S, \ldots, T \}
```

Here are Stata commands

__Step.1 Genereting Y dot__
```
xtset id year
	bysort id: gen y_dot = y - (L1.y + L2.y + L3.y)/3 if f04
	bysort id: replace y_dot = y - (L2.y + L3.y + L4.y)/3 if f05
	bysort id: replace y_dot = y - (L3.y + L4.y + L5.y)/3 if f06
```

__Step.2 Applying standard TE methods you want__

In Step 2 of Procedure 3.1, you can use built-in commands in Stata.

For example, to get Rolling RA estimates for each post-treatment period,  $t = 4, 5, 6 $, 
```
	teffects ra (y_dot x1 x2) (d) if f04, atet
	teffects ra (y_dot x1 x2) (d) if f05, atet
	teffects ra (y_dot x1 x2) (d) if f06, atet
```
For Rolling IPWRA estimates for each post-treatment period,  $t = 4, 5, 6 $
```
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f04, atet
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f05, atet
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f06, atet

```

For more details, please refer to our dofile _"lee_wooldridge_rolling_common.do"_, especially _[2] Estimation_ part


## 2. Staggered Intervension Case

__Basic set up__

1) Time Periods: $t \in {1,...,T}$

2) Time dummies: $f2_t, ...,fT_t$ , i.e., $f4 =1$ if $t=4$

3) The first Intervention occurs at $S$, $1 < S \leq T$

4) There are two observable covariates $X=(X_1, X_2)$

### Procedure 4.1
__Step.1__ For $t \in \{ g, g+1, \ldots, T \}$ and each unit $i$, compute
```math
   \dot{Y}_{igt} \equiv Y_{it}-\frac{1}{g-1} \sum_{q=1}^{g-1} Y_{iq}
```

__Step.2__ Choose as the control group the not-yet-treated units with 
```math
A_{t+1} \equiv  D_{i, t+1} + D_{i,t+2} + \cdots + D_T + D_{\infty} = 1
```

__Step.3__ Using the subset of data units ($A_{t+1} +D_g = 1$), apply standard TE methods - such as linear RA, IPW, IPWRA, PS matching - to
```math
\{ ( \dot{Y}_{igt}, D_{ig}, \mathbf{X}_i ) \quad , { i = 1, \ldots, N; g = S, \ldots, T; t = g, g+1, ..., T } \}
```

Here are Stata commands

__Step.1 Genereting Y dot__
```
xtset id year
	bysort id: gen y_dot = y - (L1.y + L2.y + L3.y)/3 if f04
	bysort id: replace y_dot = y - (L2.y + L3.y + L4.y)/3 if f05
	bysort id: replace y_dot = y - (L3.y + L4.y + L5.y)/3 if f06
```

__Step.2 Applying standard TE methods you want__

In Step 2 of Procedure 3.1, you can use built-in commands in Stata.

For example, to get Rolling RA estimates for each post-treatment period,  $t = 4, 5, 6 $, 
```
	teffects ra (y_dot x1 x2) (d) if f04, atet
	teffects ra (y_dot x1 x2) (d) if f05, atet
	teffects ra (y_dot x1 x2) (d) if f06, atet
```
For Rolling IPWRA estimates for each post-treatment period,  $t = 4, 5, 6 $
```
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f04, atet
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f05, atet
	teffects ipwra (y_dot x1 x2) (d x1 x2) if f06, atet

```
