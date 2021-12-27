%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB code developed by Helia Farzaneh (helia.farzaneh90@gmail.com) and Mohammad Karamouz (karamouz@ut.ac.ir),
% University of Tehran
%
% Last modified on September 20, 2019
%
% Please contact Helia Farzaneh with any issue.
%
% Disclaimer:
% This program (hereafter, software) is designed for instructional, educational and research use only.
% Commercial use is prohibited. The software is provided 'as is' without warranty
% of any kind, either express or implied. The software could include technical or other mistakes,
% inaccuracies or typographical errors. The use of the software is done at your own discretion and
% risk and with agreement that you will be solely responsible for any damage and that the authors
% and their affiliate institutions accept no responsibility for errors or omissions in the software
% or documentation. In no event shall the authors or their affiliate institutions be liable to you or
% any third parties for any special, indirect or consequential damages of any kind, or any damages whatsoever.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
Copulatype = 'Hic';
data=xlsread('Data1NoSandy');
[number_of_data, n] = size(data);

alpha=0.05; %p-value less than 0.05
D=data(:,1); 
L_Window=10; % Initial guess on the number of years within a window
Condition=2;
si=1;

for o=1:(length(D)-L_Window)
    Window_Data= D(1:L_Window+o-1,:);
    [H,Z,p_value]= Mann_Kendall(Window_Data(1:length(Window_Data),si),alpha);

    if H==1
        L_Window_Updated=o+L_Window-1; %Updating the number of years within one window
        Error_Counter=1;
        Window_Number=1;
        Numerator=1;

        for Numerator=1:(length(D)-L_Window_Updated+Error_Counter) %Number of windows
            Window_Data=D(Window_Number:Window_Number+L_Window_Updated-1,:);
            [H1,Z1,p_value1]= Mann_Kendall(Window_Data(1:length(Window_Data),si),alpha);

          if H1==1 %Nonstationary data in window
              Condition=1;
              L_Window_Updated=L_Window_Updated-1; %remove last data in the window to investigate the nonstationarity
              Error_function(Error_Counter,:)=[Numerator abs(p_value1)];
              Error_Counter=Error_Counter+1; %keep the numerator constant

          else % If Still stationary
              Condition=0;
              Window_Number=Window_Number+1; %move the window
          end

        end

        while Numerator<length(D)-L_Window_Updated+(Error_Counter-1) %starting from last year of previous loop

                m_Numerator=Numerator;

               for Numerator=m_Numerator:(length(D)-L_Window_Updated+(Error_Counter-1))
                   
                Window_Data=D(Window_Number:Window_Number+L_Window_Updated-1,:);
                [H1,Z1,p_value1]= Mann_Kendall(Window_Data(1:length(Window_Data),si),alpha);
                  if H1==1
                      Condition=1; % Nonstationary
                      L_Window_Updated=L_Window_Updated-1;  %remove last data in the window to investigate the nonstationarity
                      Error_function(Error_Counter,:)=[Numerator abs(p_value1)];
                      Error_Counter=Error_Counter+1;   %keep the numerator constant
                  else 
                      Condition=0; % Still stationary
                      Window_Number=Window_Number+1; %remove first data in the window
                  end
               end
        end

       if Condition==0 % Still stationary
           break
       end
    end % if H == 1

    if Condition==0 % Still stationary
        break
    end %initial for loop
end

%% L-Moment Method

L_Window=round(L_Window_Updated/10)*10;

StormSurgeStationary_gevdata=gevfit(data(:,1));
RainfallStationary_gevdata=gevfit(data(:,2));


%% Empirical copula
for m=1:number_of_data-L_Window+2
    if m<number_of_data-L_Window+2
        V=data(m:L_Window+m-1,:);
        L=V(:,1); %Surge at each window
        Surge= L;
        R=V(:,2); %rain at each window
        Precipitation= R;
    end
    if m==number_of_data-L_Window+2
        V=data;
        L=V(:,1); %Surge at each window
        Surge= L;
        R=V(:,2); %rain at each window
        Precipitation= R;
    end
        
x= [L, R]; %this is the unsorted surge rain matrix at each window
[m_Numerator, n] = size(x);

y = sort(x); %this is the sorted surge rain matrix at each window

    for i=1:m_Numerator
        for j=1:m_Numerator
            ecopula(i,j) = sum( (x(:,1)<=y(i,1)).*(x(:,2)<=y(j,2)) )/m_Numerator; %empirical copula - a matrix of 30x30
        end
    end

    for k=1:m_Numerator
     for i=1:m_Numerator
        if y(i,1) == x( k,1)
            for j=1:m_Numerator
                if y(j,2)== x(k,2)
                    H2(k,1)=ecopula(j,i); %this function indicates the copula values when sorted and unsorated data matrixes are equal which happens once per year
                    jim=[i j]; %jim determines when the surge (i) and precipitaion (j) coincide
                end
            end
        end
      end
    end

[Kendall] = corr(L,R,'type','Kendall'); % Compute Kendal's Correlation Coefficient

PACLa(m)=(2*Kendall)/(1-Kendall);
PACLa(m)=abs(PACLa(m)); %Copula parameter - Clayton
PAGUm(m)=(1)/(1-Kendall);
PAGUm(m)=max(1,PAGUm(m)); %Copula parameter - Gumbel
skipNum = 0;    
    for i= -5:0.01:5
        if ~ismember(i,skipNum)
            zart=@(t) -log(((exp(-i*t)-1))./(exp(-i)-1))./((i)./(1-exp(i*t))); %generator function in frank copula
            PF=4*(integral(zart,0,1))+1;
            oh=Kendall-PF;
            if oh>=-0.01
                if oh<=0.01
                    PAFra(m)=i; %Copula parameter - Frank
                    break
                end

            end
        end  
    end
        
    parmhat01(m,:)=gevfit(L); %parmhat(1) is the shape parameter, k, parmhat(2) is the scale parameter, sigma, and parmhat(3) is the location parameter, mu.
    [fl]=min(0.99999,gevcdf(L,(StormSurgeStationary_gevdata(:,1)),StormSurgeStationary_gevdata(:,2),parmhat01(m,3))); %fitting the location parameter in surge data - considering the location parameter be nonstationary
    parmhat02(m,:)=gevfit(R);
    [fr]=min(0.99999,gevcdf(R,(RainfallStationary_gevdata(:,1)),RainfallStationary_gevdata(:,2),RainfallStationary_gevdata(:,3))); %fitting the location parameter in rainfall data
    u=[fl fr]; % cdf matrix at each window

%% Copula Function Type        
        paramfrank(m) = copulafit('Frank',u,'Method', 'ML'); %Frank Copula - computes 99% confidence intervals for the estimated copula parameter and uses an approximation method to fit the copula
        copfrankcdf= copulacdf('Frank',u,PAFra(m)); %returns the cumulative probability of the Gaussian copula, with linear correlation parameters rho 

        paramGumbel(m) = copulafit('Gumbel',u,'Method', 'ML'); %Gumbel Copula - computes 99% confidence intervals for the estimated copula parameter and uses an approximation method to fit the copula
        copGumbelcdf= copulacdf('Gumbel',u,max(PAGUm(m),1)); %returns the cumulative probability of the Gaussian copula, with linear correlation parameters rho 

        paramClayton(m) = copulafit('Clayton',u,'Method', 'ML'); %Clayton Copula - computes 99% confidence intervals for the estimated copula parameter and uses an approximation method to fit the copula
        copClaytoncdf= copulacdf('Clayton',u,abs(PACLa(m))); %returns the cumulative probability of the Gaussian copula, with linear correlation parameters rho 

%% One parameter copulas - definition
        syms Pamh  
        amh1=((3*Pamh-1)/3*Pamh)-((2*((1-Pamh).^2)*log(1-Pamh))/3*(Pamh.^2)); % 'Ali-Mikhail-Haq' % Copula
        PAAMh(m)=double(vpasolve(amh1==Kendall,Pamh)); %PAAMh is the value where Kendall is set to equal to function above
       
        for i=1:m_Numerator
            copAMHcdf(i,:)=(u(i,1)*u(i,2))/(1-PAAMh(m)*(1-u(i,1))*(1-u(i,2))); % Based on Wikipedia: https://en.wikipedia.org/wiki/Copula_(probability_theory)
        end

        syms x1 y1 a1
        cho=(x1.*(((2*a1*x1-2*a1)*y1)-a1*x1+a1+1))*(y1.*(((2*a1*y1-2*a1)*x1)-a1*y1+a1+1)); % 'Farlie-Gumbel-Morgenstern' % Copula
        FGm1 = 4*(0.5-(int(int(cho,x1,0,1),y1,0,1)))-1;
        PFGm1(m)=double(solve(FGm1==Kendall,a1)); %PFGm1 is the value where Kendall is set to equal to function above
        for i=1:m_Numerator
            copFGMcdf(i,:)=(u(i,1)*u(i,2))*(1*PAAMh(m)*(1-u(i,1))*(1-u(i,2))); % Based on: https://people.math.ethz.ch/~embrecht/ftp/copchapter.pdf & http://www.maths.manchester.ac.uk/~saralees/chap20.pdf
        end

    PARAman(m,:)=[Kendall PACLa(m) paramClayton(m) PAGUm(m) paramGumbel(m) PAFra(m) paramfrank(m) PAAMh(m) PFGm1(m)];

    % tail is equal to 0.05 and data size is 30
    KSthreshold = 1.358/(m_Numerator)^0.5;
    
    % testing each copula function using KS threshold
    KSfrankvalue(m) = max (sort(H2)-sort (copfrankcdf)); 
    if  KSfrankvalue(m) < KSthreshold
       OLSfrankvalue(m) =(sum ((H2-copfrankcdf).^2)/m_Numerator).^0.5;
    else 
     OLSfrankvalue(m)=10;
    end 
    
    KSGumbelvalue(m)= max (sort(H2)-sort (copGumbelcdf));
    if KSGumbelvalue(m) < KSthreshold
       OLSGumbelvalue(m)=(sum ((H2-copGumbelcdf).^2)/m_Numerator).^0.5;
       else
       OLSGumbelvalue(m)=10;
    end

    KSClaytonvalue(m)= max (sort(H2)-sort (copClaytoncdf));
    if KSClaytonvalue(m) < KSthreshold
       OLSClaytonvalue(m)=(sum ((H2-copClaytoncdf).^2)/m_Numerator).^0.5;
       else
       OLSClaytonvalue(m)=10;
    end
    
    KSAMHvalue(m)= max (sort(H2)-sort (copAMHcdf));
    if KSAMHvalue(m) < KSthreshold
       OLSAMHvalue(m)=(sum ((H2-copAMHcdf).^2)/m_Numerator).^0.5;
    else
        OLSAMHvalue(m)=10;
    end
    
    KSFGMvalue(m)= max (sort(H2)-sort (copFGMcdf));
    if KSFGMvalue(m) < KSthreshold
       OLSFGMvalue(m)=(sum ((H2-copFGMcdf).^2)/m_Numerator).^0.5;
    else
       OLSFGMvalue(m)=10;
    end
     
    % testing each copula function using OLS value
    OLS(m,:) = [ OLSClaytonvalue(m)+10, OLSGumbelvalue(m), OLSfrankvalue(m)];
    [M,I]= min (OLS(m,:));
   
    switch I
        case 1
            fprintf('best copula is Clayton for window %d\n',m)
            PARfinal(m,:)=[PACLa(m) 1];
            id='clayton';
        case 2
            fprintf('best copula is Gumbel for window %d\n',m )
            PARfinal(m,:)=[PAGUm(m) 2];
             id='gumbel';
        otherwise
            fprintf('best copula is frank for window %d\n',m)
            PARfinal(m,:)=[PAFra(m) 3];
            id='frank';

    end

%% Figure
meshspacing = [linspace(1e-5,.95,1000) linspace(0.95+1e-5,1-1e-5,2000)];
[xx,yy] = meshgrid(meshspacing,meshspacing); %Defining the plotting space for each figure

ToleranceL=10*(max(L)-min(L))/m_Numerator; %Margin for surge is defiend to extend the lower and upper limits
S1=linspace(min(L)-ToleranceL,max(L)+ToleranceL,3000);
[xxs1,yys1] = meshgrid(S1,S1); %creating mesh grid for data point generation

parmhat511(m,:)=gevfit(L); %gev parameters for each time window of 30 years - surge

Surge_CDF_Each_Window=gevcdf(L,(StormSurgeStationary_gevdata(:,1)),StormSurgeStationary_gevdata(:,2), parmhat511(m,3)); %fitting the gev parameters to the surge data for calculating the cdf function - nonstationarity is considered in location parameter of surge only
Precipitation_CDF_Each_Window=gevcdf(R,(RainfallStationary_gevdata(:,1)),RainfallStationary_gevdata(:,2),RainfallStationary_gevdata(:,3)); %fitting the gev parameters to the rainfall data for calculating the cdf function
CDF_Matrix_Each_Window=[Surge_CDF_Each_Window Precipitation_CDF_Each_Window];
Copula_CDF_Each_Window=copulacdf(id,CDF_Matrix_Each_Window, PARfinal(m,1)); % Copula cumulative distribution function
C00=1-Surge_CDF_Each_Window-Precipitation_CDF_Each_Window+Copula_CDF_Each_Window; %1-cdf of surge-cdf of rainfall +copulacdf
Return_Interval_Each_Year=1./(C00); %Return Intervals for 30 years within each window 30x1 matrix

A=min(0.999,gevcdf(xxs1,(StormSurgeStationary_gevdata(:,1)),StormSurgeStationary_gevdata(:,2),parmhat511(m,3))); %Now plotting the points on the figure - remember A points are for surge
A=reshape(A,[(size(xxs1,1)*size(xxs1,1)),1]); %reshaping A into column matrix

% We now repeat every single step for precipiation 
ToleranceR=10*(max(R)-min(R))/m_Numerator; %Margin for precipitation is defiend to extend the lower and upper limits
S2=linspace(min(R)-ToleranceR,max(R)+ToleranceR,3000);
[xxs2,yys2] = meshgrid(S2,S2); %creating mesh grid for data point generation

B=min(0.999,gevcdf(yys2,(RainfallStationary_gevdata(:,1)),RainfallStationary_gevdata(:,2),RainfallStationary_gevdata(:,3))); %Now plotting the points on the figure - remember B points are for rainfall
B=reshape(B,[(size(xxs1,1)*size(xxs1,1)),1]); %reshaping A into column matrix

S = [A,B]; %Matrix of new surge and precipitation data 2 columns
copidCDF=copulacdf(id,S, PARfinal(m,1)); % Copula cumulative distribution function

C=1-A-B+copidCDF;
T=1./((C));
P = [25]; %25-years is selected as the return period
Tolerance=0.005;
P_LB = P - Tolerance*P; P_UB = P + Tolerance*P;
[P_Sort, ID_Pr] = sort(T);

Font_SIZE=10;

for j = 1:1
    % Find indices associated with each probability contour
    ID_Contour = find( P_Sort(:) >= P_LB(j) & P_Sort(:) <= P_UB(j) );
    Countors{j}=ID_Contour;
   
    
    % Sort first uniform marginal
    U1  = S(ID_Pr, 1);
    U2  = S(ID_Pr, 2);
    [UU, ID_U] = sort( U1(ID_Contour) );
    
    % Associated second uniform marginal
    VV = U2(ID_Contour);
    VVV = VV(ID_U);

    %Estimate the actual variable value associated with the probability
    parmhat31(m,:)=gevfit(L);
    IUU=gevinv(UU,( StormSurgeStationary_gevdata(:,1)), StormSurgeStationary_gevdata(:,2),parmhat31(m,3));
 
    parmhat32(m,:)=gevfit(R);
    IVVV=gevinv(VVV,(RainfallStationary_gevdata(:,1)),RainfallStationary_gevdata(:,2),RainfallStationary_gevdata(:,3));
    Points{j}=[IUU IVVV];  %data points to plot
       
    % Trick: remove infinity if there is any
    IDNAN = find( isinf(IUU) | isinf(IVVV) );
    IUU(IDNAN) = []; IVVV(IDNAN) = [];
    
    % Plot each probability contour
    %if length(IUU) > 4
    hold on
    if m<number_of_data-L_Window+2
    [maxPDFvalue(m,1),valueS(m,1)]=max(copulapdf(id,[UU VVV],PARfinal(m,1)));
    designValueSurge(m,:)=gevinv(UU(valueS(m,1),1),( StormSurgeStationary_gevdata(:,1)), StormSurgeStationary_gevdata(:,2), parmhat31(m,3));
    designValuePrecipitation(m,:)=gevinv(VVV(valueS(m,1),1),(RainfallStationary_gevdata(:,1)),RainfallStationary_gevdata(:,2),RainfallStationary_gevdata(:,3));
        
    Co1(m,:)=[((m-1)/(number_of_data-L_Window)) 0 1-(((m-1)/(number_of_data-L_Window)))];
    TemporalCopulaValue=plot(IUU, IVVV, 'color', Co1(m,:), 'linewidth',0.5);
    Year(m,:) = num2str(2017-(number_of_data-L_Window)+m-1);
    end
    
    if m==number_of_data-L_Window+2
        [maxPDFvalue(m,1),valueS(m,1)]=max(copulapdf(id,[UU VVV],PARfinal(m,1)));
        designValueSurge(m,:)=gevinv(UU(valueS(m,1),1),( StormSurgeStationary_gevdata(:,1)), StormSurgeStationary_gevdata(:,2),StormSurgeStationary_gevdata(:,3));
        designValuePrecipitation(m,:)=gevinv(VVV(valueS(m,1),1),(RainfallStationary_gevdata(:,1)),RainfallStationary_gevdata(:,2),RainfallStationary_gevdata(:,3));
        Co1(m,:)=[0 0 0];
        TemporalCopulaValue=plot(IUU, IVVV, 'color', Co1(m,:), 'linewidth',10);
        Year(m,:) = 'Stat'
    end
    
    [hleg1, hobj1] = legend(Year);
    textobj = findobj(hobj1, 'type', 'text');
    set(textobj, 'Interpreter', 'latex', 'fontsize', Font_SIZE);
%     legend('Location','eastoutside')
    
%     hleg1.Location = 'eastoutside';
%     HeightScaleFactor = 1.5;
%     NewHeight = hleg1.Position(4) / HeightScaleFactor;
%     hleg1.Position(2) = hleg1.Position(2) - (NewHeight - hleg1.Position(4));
%     hleg1.Position(4) = NewHeight;
    
    hold on
    xlabel ('Storm surge (m)')
    ylabel ('Precipitation (mm)')
    
    t = title(strcat('\color{red}',upper(id)));
    set(t, 'horizontalAlignment', 'right','fontname','times','fontweight','bold','fontsize',Font_SIZE+4);
    set(t, 'units', 'normalized');
    h1 = get(t, 'position');
    set(t, 'position', [1 h1(2) h1(3)])
    hold on
    drawnow
end
end


