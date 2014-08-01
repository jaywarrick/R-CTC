R1 <- seq(0,1,0.001)
R2 <- 1
PI <- pi

getMin <- function(R1, R2)
{
  return(((4*pi*(R2^3/3-(R2^2-R1^2)^(3/2)/3)-(4/3)*pi*R1^3)/((4/3)*pi*R2^3-(4/3)*pi*R1^3)))
}

plot(R1,getMin(R1,R2), type='l')

Cyt <- 0.2  #Conc=0.5*(((4/3)*pi*R2^3) - ((4/3)*pi*R1^3))
Nuc <- 0.8

ObsNucSig <- Nuc + getMin(R1,R2)*Cyt
ObsCytSig <- (1-getMin(R1,R2))*Cyt

Tot <- ObsNucSig + ObsCytSig

Ratio <- ObsNucSig/Tot

plot(R1,Ratio)

plot(R1,getMin(R1,R2), type='l')

MeasuredNucPercentage <- .7 # (i.e., Nuc Sig / Total Signal)
R1 <- 0.8
R2 <- 1
ActualNucPercentage <- (MeasuredNucPercentage - getMin(R1, R2))/(1-getMin(R1,R2))