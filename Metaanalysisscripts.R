### Script for inverse variance weighted meta-analysis
## Inputs: Bi and SEi

weight = function(SE){
  weight_i = 1/(SE^2)
  return(weight_i)
}


### Meta-analysis of LEAP and CHOP

BetaCHOP = 0.036
SECHOP = 0.13
WeightCHOP = weight(SECHOP)

BetaLEAP = 0.0199
SELEAP = 0.051
WeightLEAP = weight(SELEAP)

SEtotal = sqrt(1/(WeightLEAP + WeightCHOP))
Betatotal = (BetaLEAP*WeightLEAP + BetaCHOP*WeightCHOP)/(WeightLEAP + WeightCHOP)
Ztotal = Betatotal/SEtotal

Print(Betatotal)
Print(SEtotal)
Print(Ztotal)

#Betatotal = 0.022
#SEtotal = 0.047
#Ztotal = 0.46

### Meta-analysis of LEAP, CHOP and SSC

BetaSSC =  0.052
SESSC = 0.020
WeightSSC = weight(SESSC)



SEtotal = sqrt(1/(WeightLEAP + WeightCHOP + WeightSSC))
Betatotal = (BetaLEAP*WeightLEAP + BetaCHOP*WeightCHOP + BetaSSC*WeightSSC)/(WeightLEAP + WeightCHOP + WeightSSC)
Ztotal = Betatotal/SEtotal

print(Betatotal) # 0.04748581
print(SEtotal) # 0.01843137
print(Ztotal) # 2.576357

### ADOS meta in CHOP and SSC
BetaCHOP = -0.046
SECHOP = 0.04
WeightCHOP = weight(SECHOP)

BetaSSC = -0.00099
SESSC = 0.018
WeightSSC = weight(SESSC)

SEtotal = sqrt(1/(WeightSSC + WeightCHOP))
Betatotal = (BetaSSC*WeightSSC + BetaCHOP*WeightCHOP)/(WeightSSC + WeightCHOP)
Ztotal = Betatotal/SEtotal


# Betatotal = -0.008569647
# SEtotal = 0.01641459
# Ztotal = -0.5220751







