#R version 3.6.1, not Rstudio
rm(list=ls())
library(nleqslv) #library for solving nonlinear equations

#³æ¦ì million dallar(indonesia's equities need to be divided by 1000)
#obtain data from the file in the form of .csv by the path
#using processed data
#data=as.matrix(read.csv("C:/Users/user/Documents/IN_ID_HK_TW_CN2020-01-01(2)processed.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/IN_ID_HK_TW_CN20191231processed.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/IN_ID_HK_TW_CN20191001processed.csv",head=FALSE))
#data=as.matrix(read.csv("C:/Users/user/Documents/IN_ID_HK_TW_CN20190701processedup.csv",head=FALSE))#6mon
#==========================================================================

#using 2 sets of raw data
#raw_data1=as.matrix(read.csv("C:/Users/user/Documents/rawdata1.csv",head=FALSE))
#raw_data2=as.matrix(read.csv("C:/Users/user/Documents/rawdata2.csv",head=FALSE))
raw_data1=as.matrix(read.csv("C:/Users/user/Documents/rawdata_new_down1.csv",head=FALSE))
raw_data2=as.matrix(read.csv("C:/Users/user/Documents/rawdata_new_down2.csv",head=FALSE))
data=matrix(0,NROW(raw_data1)-4,3)#save result

for(i in 1: (NROW(data)/3) ){
  num1=3*i-2; num2=3*i-1; num3=3*i
  #P1
  E_0=as.numeric(raw_data1[num1+4,4])
  E_t=as.numeric(raw_data2[num1+4,4])
  sigma_e=sd(raw_data1[num2+4,c(4:NCOL(raw_data1))],na.rm = TRUE)/100*(260^0.5)
  #P2
  sigma_et=sd(raw_data2[num2+4,c(4:NCOL(raw_data2))],na.rm = TRUE)/100*(260^0.5)
  L_t=as.numeric(raw_data1[num3+4,NCOL(raw_data1)])
  L_t2=as.numeric(raw_data2[num3+4,NCOL(raw_data2)])
  
  # when number>10^8, numerical method usually cannot find solutions.
  # use unit: million
  if(!(is.na(E_0)|is.na(E_t))){
    if(E_0>10^8|E_t>10^8){
      data[num1,2]=E_0/1000
      data[num1,3]=E_t/1000
      data[num3,2]=L_t/1000
      data[num3,3]=L_t2/1000
    }
    else{
      data[num1,2]=E_0
      data[num1,3]=E_t
      data[num3,2]=L_t
      data[num3,3]=L_t2
    }
  }
  else{
    data[num1,2]=E_0
    data[num1,3]=E_t
    data[num3,2]=L_t
    data[num3,3]=L_t2
  }
  
  data[num2,2]=sigma_e
  data[num2,3]=sigma_et
}
#name
data[ ,1]=raw_data1[5:NROW(raw_data1) ,1]

#==========================================================================

#end-function: solve two nonlinear equations for v_a; sigma_a 
#newton numerical method
#the start date is regared as time point 0; the end date as time point t.
#therefor, the duration equals t-0=t
#E_0 means the total equity market value at time point 0
#sigma_e means the standar deviation of the equity market value of one year
#L_t means total liebilities at time piont t
#risk_f means risk free rate
#t= duration; name= name of the company
end<- function( E_0, sigma_e, L_t,  risk_f, t, name, sigma_e0=0){
  
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
  p_edf=0
  if (sigma_e0!=0){
    #adjust DD: (*sigma_et/sigma_e0)!!!!
    DD=(log(v_0)-log(L_t)+(mu_v-1/2*sigma_v^2)*t)/(sigma_v*sqrt(t))*sigma_et/sigma_e0
    p_edf=pnorm(-DD,0,1)
  }
  else{
    DD=(log(v_0)-log(L_t)+(mu_v-1/2*sigma_v^2)*t)/(sigma_v*sqrt(t))
    p_edf=pnorm(-DD,0,1)
  }

  
  #warnning level
  default_risk=FALSE
  if (p_edf>=0.01){
    default_risk= TRUE
    print(c("name=",name,"v_0=",v_0,"sigma_vt=",sigma_v,"p_edf",p_edf))
  }
  
  name=name
  return(matrix(c(name,-1*DD,result[["message"]],default_risk)))
  
}

#==================================================================================
#create matrix to save returns of function called
results_all=matrix(0,NROW(data)/3,9)

#input data to function that is named end
#end( E_0, sigma_e, L_t,  risk_f, t, name, adjust sigma)
for(i in 1:(NROW(data)/3)){
  num1=3*i-2;num2=3*i-1;num3=3*i
  name=as.character(data[num1,1])
  
  if (!(name %in% results_all[,1])){  # if the companies are repeat, do not run again
    #p1
    E_0=as.numeric(data[num1,2]); sigma_e=(3/3)*as.numeric(data[num2,2])+(0/3)*as.numeric(data[num2,3])
    L_t=as.numeric(data[num3,2]); 
    if(is.na(E_0)| is.na(sigma_e)| is.na(L_t)){
      results_all[i,1:2]=c(name,'data missing')
    }
    else{
      #first period need not be adjusted(pass 0 to function)
      results_all[i,1:4]=end(E_0,sigma_e,L_t,0.01,3/12,name,0) 
    }
    
    #p2
    E_t=as.numeric(data[num1,3]); sigma_et=(0/3)*as.numeric(data[num2,2])+(3/3)*as.numeric(data[num2,3])
    L_t2=as.numeric(data[num3,3])
    if(is.na(E_t)| is.na(sigma_et)| is.na(L_t2)){
      results_all[i,5:6]=c(name,'data missing')
    }
    else{
      #second period need to be adjust(pass the sigma of last period: sigma_e)
      results_all[i,5:8]=end(E_t,sigma_et,L_t2,0.01,3/12,name,sigma_e)
    }
   
    #compare
    p1=as.numeric(results_all[i,2]);p2=as.numeric(results_all[i,6])
    #both P1 and P2 do not miss
    if(!(is.na(p1) | is.na(p2))){
      
      if(p1<p2){
        results_all[i,9]="down"
      }
      else if(p1>p2){
        results_all[i,9]="up"
      }
      else{
        results_all[i,9]="not change"
      }
    }
    #P1 or P2 miss
    else{
      results_all[i,9]="p1 or p2 missing"
    }
  }
}

#delete the empty rows that is caused by repeat
results_all=subset(results_all,results_all[,1]!=0)  


#===============================================================

# output as csv, filtered by some condition: result_all
# default_company= matrix(c('name','movement'),1,2)
# for(i in 1:NROW(results_all)){
#   if (results_all[i,9]=="up"){
#     temp=c(results_all[i,1],results_all[i,9])
#     #conbine the temp to the matrix by row
#     default_company=rbind(default_company,temp)
#   }
# }
# write.csv(default_company,"auto.csv",row.names = FALSE)
#==============================================================

# output as csv, filtered by some condition: data missing
company= matrix(c('name','total market value','price change 1D','total liability'),1,4)
for(i in 1:(NROW(data)/3)){
  num1=3*i-2;num2=3*i-1;num3=3*i
  
  name=as.character(data[num1,1])
  
  E_0=0
  sigma_e=0
  L_t=0
  E_0=as.numeric(data[num1,2])
  sigma_e=as.numeric(data[num2,2])
  L_t=as.numeric(data[num3,2])
  
  data_missing=c(0,0,0,0)
  if(is.na(E_0)| is.na(sigma_e)| is.na(L_t)){
    data_missing[1]=name
    if(is.na(E_0)){
      data_missing[2]='True'
    }
    if(is.na(sigma_e)){
      data_missing[3]='True'
    }
    if(is.na(L_t)){
      data_missing[4]='True'
    }
    if (!(name %in% company[,1])){
      company=rbind(company,data_missing)
    }
    
  }
  
}
write.csv(company,"auto.csv",row.names = FALSE)


