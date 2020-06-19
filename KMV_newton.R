rm(list=ls())
library(nleqslv)

#³æ¦ì million dallar
#obtain data from the file in the form of .csv by the path
#data=as.matrix(read.csv("C:/Users/user/Documents/tw900-processed.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/JCI700_ID_processed.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/IN_500_processed.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/HK_RATC_D_P.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/IN_RATC_D_P.csv",head=FALSE))
data=as.matrix(read.csv("C:/Users/user/Documents/ID_RATC_D_P.csv",head=FALSE))

#data[1,1]=2330

#solve two nonlinear equations for v_a; sigma_a 
#newton numerical method
#the start date is regared as time point 0; the end date as time point t.
#therefor, the duration equals t-0=t
#E_0 means the total equity market value at time point 0
#sigma_e means the standar deviation of the equity market value of one year
#L_t means total liebilities at time piont t
#risk_f means risk free rate
#t= duration; name= name of the company
end<- function( E_0, sigma_e, L_t,  risk_f, t, name){
  
  KMV<- function(roots){
    f<- numeric(length (roots))
    
    d1=((log(roots[1])-log(L_t)+(risk_f+0.5*roots[2]^2)*t))/(roots[2]*sqrt(t))
    d2=((log(roots[1])-log(L_t)+(risk_f-0.5*roots[2]^2)*t))/(roots[2]*sqrt(t))
    
    f[1]= roots[1]*pnorm(d1)-
      L_t*exp(-risk_f*t)*pnorm(d2)-
      E_0
    f[2]=(roots[1]*roots[2]/E_0)*pnorm(d1)-
      sigma_e
    f
  }
  
  #initial value
  #start_root<-c(v_e,sigma_e)
  start_root<-c(E_0+L_t*exp(-risk_f*t),sigma_e*(E_0/(E_0+L_t*exp(-risk_f*t))))
  
  result= nleqslv(start_root,KMV, method = "Newton", control = list(allowSingular=TRUE, maxit=2500))
  
  #calculate p-edf
  #v_0 means the real value of total assets derived from KMV
  #sigma_v means the standar deviation derived from KMV
  v_0=result$x[1];sigma_v=result$x[2]
  #mu_v means groth rate of assets within one year
  mu_v=0
  DD=(log(v_0)-log(L_t)+(mu_v-1/2*sigma_v^2)*t)/(sigma_v*sqrt(t))
  p_edf=pnorm(-DD,0,1)
  
  #warnning level
  default_risk=FALSE
  if (p_edf>=0.01){
    default_risk= TRUE
    print(c("name=",name,"v_0=",v_0,"sigma_vt=",sigma_v,"p_edf",p_edf))
  }
  
  name=name
  return(matrix(c(name,p_edf,result[["message"]],default_risk)))
  
}

#create matrix to save returns of function called
name_pedf_newton=matrix(0,NROW(data)/3,4)

#import data to function named end
#end( E_0, sigma_e, L_t,  risk_f, t, name)
for(i in 1:(NROW(data)/3)){
  num1=3*i-2;num2=3*i-1;num3=3*i
  E_0=as.numeric(data[num1,2])/1000; sigma_e=as.numeric(data[num2,2])
  L_t=as.numeric(data[num3,2])/1000; name=data[num1,1]
  if(is.na(E_0)| is.na(sigma_e)| is.na(L_t)){
    name_pedf_newton[i,1:2]=c(name,'data missing')
  }
  else{
    name_pedf_newton[i,1:4]=end(E_0,sigma_e,L_t,0.01,1,name)
  }
  
}

# default_company= matrix(c('name','probability'),1,2)
# for(i in 1:NROW(name_pedf_newton)){
#   if (name_pedf_newton[i,4]== TRUE){
#     temp=c(name_pedf_newton[i,1],name_pedf_newton[i,2])
#     default_company=rbind(default_company,temp)
#   }
# }

#write.csv(default_company,"auto.csv",row.names = FALSE)
