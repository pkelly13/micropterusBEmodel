#Bioenergetics model for largemouth bass - adapted from Kitchell and Stewart 1977 and Rice et al. 1983
#Patrick Kelly 15 July 2014

#overarching model is dBbdt = C - (Rsa + Rsda +F +U) in kJ g-1 d-1 - conversion rate -> 1g = 4.184 kJ to convert to biomass rate

#start by modeling consumption
#C = Cmax * rc * P
#Cmax = a1*B*b1

a1<-0.33 #intercept for maximum consumption (g g-1 day-1)
b1<--0.325 #weight dependence exponent for maximum consumption
B<-400 #biomass in grams - user input <--------------

Cmax<-a1*B^b1

Tm<-37 #maximum temperature (from Rice et al 1983)
Topt<-27.5 #optimal temperature (from Rice et al 1983)
Tamb<-19 #user input <-------------
V<-(Tm-Tamb)/(Tm-Topt)

Q<-2.65 #slope value from Rice et al. 1983
W<-log(Q)*(Tm-Topt)
Y<-log(Q)*(Tm-Topt+2)
X<-((W^2*(1+(1+40/Y)^0.5)^2))/400

rc<-(V^X)*exp(X*(1-V)) #temperature dependent proportional adjustment (0 to 1) of consumption rate - from Kitchell and Stewart 1977

P<-0.7 #proportionality constant from 0 to 1, used to adjust the ration when fitting an observed growth curve - user input <---------------

C<-Cmax*rc*P #consumption rate - could think of a different way to model consumption rates perhaps

#Respiration - equation from Rice et al 1983 doesn't make sense, so used the equation from 
#R = Rmax * A * rR + SC

a2<-0.348 #intercept for respiration (mg O2 g-1 hr-1)
b2<--0.355 #weight dependence exponent for repsiration
Rmax<-a2*B^b2

A=1 #activity coefficient (from Kitchell 1977)

Topt<-28
Tm<-34
Q<-2.1
W<-log(Q)*(Tm-Topt)
Y<-log(Q)*(Tm-Topt+2)
X<-((W^2*(1+(1+40/Y)^0.5)^2))/400

rR<-(V^X)*exp(X*(1-V)) #a temperature dependent proportional adjustment of repsiration

S<-0.142 #proportion of consumed energy utilixed in apparent SDA

R<-Rmax*A*rR+(S*C) #respiration/metabolism

#Egestion
f<-0.104 #proportion of consumed energy lost through egestion 
F<-f*C

#Excretion
u<-0.079 #proportion of energy lost through excretion
U<-u*C

dB.Bdt<-C-(R+F+U) #specific growth rate equation



