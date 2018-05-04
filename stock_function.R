#########################
##                     ##
##   Stock Abundance   ##
##   and Recruitment   ##
##                     ##
#########################


rho.R<-0.5
sigma.R<-0.5

## specifies constants required for calculating recruitments


stock.abundance.and.recruitment<-function(model.type, arb.val,
          sig.r, ro.r, mats, wts, MORTS, F.mortality, rho.mort, fish.age, period){

## sets up function, utilising model choice, an arbitrary value
## for age1/year1, recruitment constants sigma R and rho R,
## matrices of maturity, weight, total and fishing mortality at age/
## per year, fish ages and the years we are interested in.

source("C:/Work/Coby's Stuff/mortality functions.txt")
source("C:/Work/Coby's Stuff/weights.txt")
source("C:/Work/Coby's Stuff/Maturity.txt")
source("C:/Work/Coby's Stuff/model_selection_function.txt")

#  model.type <- "Power"
#  arb.val  <- 10000
#  sig.r    <- 0.5
#  ro.r     <- 0.5
#  mats     <- maturity.at.age.by.year
#  wts      <- weight.at.age.by.year
#  MORTS    <- TOTAL.MORTALITY
#  F.mortality <- fishing.mortality
#  rho.mort <-0.25
#  fish.age <- age
#  period   <- years
  
  
## reads in required files

stock.size<-matrix(nrow=length(period), ncol=length(fish.age))

## sets up appropriate sized matrix to recieve abundance data


sigma.d.R<-sig.r*(sqrt(1-ro.r^2))

z.y<-sigma.d.R*rnorm(length(period), mean=0, sd=1)

epsilon.y.R<-vector(length=length(period))

for (i in (1:length(period))){
  epsilon.y.R[i]<-ifelse(i==1, 0, ro.r*epsilon.y.R[i-1]+z.y[i])
}

## calculates recruitment variables which simplify things
## later on

expected.recruitment<-vector(length=length(period))
sp.st.bio<-vector(length=length(period))

expected.recruitment[1]<-arb.val
stock.size[1,1]<-arb.val

## sets up vectors for our expected recruitment and SSB values,
## and inserts our arbitrary value in the appropriate places

##############
##  Year 1  ##
##############

for (h in (2:length(fish.age))){
stock.size[1,h]<-stock.size[1,h-1]*exp(-MORTS[1,h-1])
}

## uses arbitrary number of age1/year1 fish and total mortality
## values to calculate numbers of fish in other age classes in first
## year natural mortality in first year and arbitrary value to
## estimate stock size in first year


sp.st.bio[1]<-sum(stock.size[1,]*mats[1,]*wts[1,])

## sums product of abundance, maturity ogive and weight at age
## to give spawning stock biomass value




##############
##  Year 2  ##
##############

if (model.type=="Ricker")
stock.size[2,1]<-((5.436*sp.st.bio[1])*exp(-0.0002*sp.st.bio[1]))
if (model.type=="Beverton-Holt")
stock.size[2,1]<-((15000*sp.st.bio[1])/2500+sp.st.bio[1])
if (model.type=="Power")
stock.size[2,1]<-1189.2071*sp.st.bio[1]^0.25
if (model.type=="Saila-Lorda")
stock.size[2,1]<-(0.002956*sp.st.bio[1]^2)*exp(-0.0004*sp.st.bio[1])
if (model.type=="Shepherd")
stock.size[2,1]<-(4*sp.st.bio[1])/(1+((sp.st.bio[1]/5000)^2))
if (model.type=="Mixed")
stock.size[2,1]<-0.001772*(sp.st.bio[1]^2)*exp(-0.0004646*sp.st.bio[1])+
            ((6665.676*sp.st.bio[1])/(2152.310+sp.st.bio[1]))+1000
if (model.type=="Changepoint")
stock.size[2,1]<-ifelse(sp.st.bio[1]<5000, ((10000*sp.st.bio[1])/5000), 10000)

expected.recruitment[2]<-stock.size[2,1]

## takes SSB value for first year and applies the chosen model
## in this case, actual recruitment == expected recruitment

stock.size[2,2]<-stock.size[1,1]*exp(-(MORTS[1,1]-F.mortality[1,1]))

for (t in (2:(length(fish.age)-1))){
stock.size[2,t]<-stock.size[1,t-1]*exp(-MORTS[1,t-1])
}

stock.size[2,length(fish.age)]<-stock.size[1,(length(fish.age)-1)]*exp(-MORTS[1,(length(fish.age)-1)])+
                            stock.size[1,(length(fish.age))]*exp(-MORTS[1,(length(fish.age))])

## generates abundance values for other year classes based on the
## correct cohort and natural mortality values

sp.st.bio[2]<-sum(stock.size[2,]*mats[2,]*wts[2,])

## produces a SSB value for our second year



##################
##  Years 3->n  ##
##################

for (i in (3:length(period))){

  if (model.type=="Ricker")
    
    expected.recruitment[i]<-((5.436*sp.st.bio[i-1])*exp(-0.0002*sp.st.bio[i-1]))
  
  if (model.type=="Beverton-Holt")
    
    expected.recruitment[i]<-((15000*sp.st.bio[i-1])/(2500+sp.st.bio[i-1]))
  
  if (model.type=="Power")
    
    expected.recruitment[i]<-1189.2071*sp.st.bio[i-1]^0.25
  
  if (model.type=="Saila-Lorda")
    
    expected.recruitment[i]<-(0.002956*sp.st.bio[i-1]^2)*exp(-0.0004*sp.st.bio[i-1])
  
  if (model.type=="Shepherd")
    
    expected.recruitment[i]<-(4*sp.st.bio[i-1])/(1+((sp.st.bio[i-1]/5000)^2))
  
  if (model.type=="Mixed")
    
    expected.recruitment[i]<-0.001772*(sp.st.bio[i-1]^2)*exp(-0.0004646*sp.st.bio[i-1])+
      ((6665.676*sp.st.bio[i-1])/(2152.310+sp.st.bio[i-1]))+1000
  
  if (model.type=="Changepoint")
    
    expected.recruitment[i]<-ifelse(sp.st.bio[i-1]<5000, ((10000*sp.st.bio[i-1])/5000), 10000)
  
  
    stock.size[i,1]<-expected.recruitment[i]*exp(epsilon.y.R[i]-(0.5*sig.r^2))
    stock.size[i,2]<-stock.size[i-1,1]*exp(-(MORTS[i-1,1]-F.mortality[i-1,1]))
    stock.size[i,length(fish.age)]<-stock.size[i-1,length(fish.age)-1]*exp(-MORTS[i-1,length(fish.age)-1])+
    stock.size[i-1,length(fish.age)]*exp(-MORTS[i-1,length(fish.age)])    

    for (j in (3:(length(fish.age)-1))){


    stock.size[i,j]<-stock.size[i-1,j-1]*exp(-MORTS[i-1,j-1])


}

## Takes SSB for the previous year and calculates expected recruitment
## based on the model choice, before calculating actual recruitment,
## then applying natural mortalities to cohorts to arrive at abundances
## for other year classes.


sp.st.bio[i]<-sum(stock.size[i,]*mats[i,]*wts[i,])

## calculates SSB for each year


}

rownames(stock.size)<-period
colnames(stock.size)<-fish.age

## makes our abundance matrix look nice

catch.numbers<-matrix(nrow=length(period), ncol=length(fish.age))
catch.numbers[,1]<-0

for (e in (1:length(period))){
    for (f in (2:length(fish.age))){

    catch.numbers[e,f]<-(stock.size[e,f]*F.mortality[e,f]*(1-exp(-MORTS[e,f])))/
        MORTS[e,f]
    }
}


total.biomass<-vector(length=length(period))
mean.F<-vector(length=length(period))
true.yield<-vector(length=length(period))
reported.yield<-vector(length=length(period))
mean.total.mortality<-vector(length=length(period))

for (p in (1:length(period))){
total.biomass[p]<-sum(stock.size[p,]*wts[p,])
mean.F[p]<-mean(F.mortality[p,2:4])
mean.total.mortality[p]<-mean(MORTS[p,2:4])
true.yield[p]<-sum(catch.numbers[p,]*wts[p,])
reported.yield[p]<-ifelse(period[p]<2015, true.yield[p], (1-rho.mort)*true.yield[p])
}

rownames(catch.numbers)<-rownames(stock.size)
colnames(catch.numbers)<-colnames(stock.size)

stock.info<-list(stock.size=stock.size, model=model.type, years=period, expected.recruitment=expected.recruitment,
actual.recruitment=as.vector(stock.size[,1]), SSB=sp.st.bio, total.biomass=total.biomass,
mean.F=mean.F, mean.total.mortality=mean.total.mortality, true.yield=true.yield, landings=round(catch.numbers), reported.yield=reported.yield,
F.mortality=F.mortality, lengths=length.at.age.by.year, weights=weight.at.age.by.year)

## combines stock size, model choice, expected and realised recruitments into a list

return(stock.info)

## returns all the stuff we're interested in

}

sim.stock<-stock.abundance.and.recruitment(model, 10000,
          sigma.R, rho.R, maturity.at.age.by.year, weight.at.age.by.year, TOTAL.MORTALITY,
          fishing.mortality, 0.25, age, years)

sim.stock
