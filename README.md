## Lee and Wooldridge 2023
_"A Simple Transformation Approach to Difference-in-Differences Estimation for Panel Data"_  
	Available on SSRN https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4516518

## How to Apply Rolling Method to the Panal Data
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
\{ ( \dot{Y}_{it}, D_i, \mathbf{X}_i) \ : \ i \ = \ 1, \ldots, N , \quad t= S, \ldots, T \}
```

Here are Stata commands

__Step.1 Genereting Y dot__
```

```
