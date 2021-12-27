function [H,Z,p_value]= Mann_Kendall(V,alpha)
V= reshape(V,length(V),1); 
alpha= alpha/2; 
n= length(V); 
i=0; j=0; S=0; 
for i=1:n-1
   for j= i+1:n 
      S= S + sign(V(j)-V(i)); 
   end
end
VarS= (n*(n-1)*(2*n+5))/18;
StdS= sqrt(VarS); 
if S >= 0
   Z= ((S-1)/StdS)*(S~=0);
else
   Z= (S+1)/StdS;
end
p_value= 2*(1-normcdf(abs(Z),0,1)); 
pz= norminv(1-alpha,0,1); 
H= abs(Z)>pz;  
return
