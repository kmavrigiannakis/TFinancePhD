%function [betaCoeff, betaStError] = betaEstimation( dependentVar, independentVar )
clear all;
load( 'monthlyReturns.mat' );
    %Remove data for 2018
    monthlyReturns = monthlyReturns1(1:size(monthlyReturns1,1)-1 ,:);
    %Correct data for NaN elements with value -99.99
    monthlyReturns( monthlyReturns==-99.99 ) = NaN;
    %monthlyReturns=monthlyReturns/100;

load('marketMinusRF.mat');
load('RF.mat');

marketMinusRF = marketMinusRF(1:1086,1);
excessReturn = monthlyReturns - RF;

independentVar = reshape( marketMinusRF, 1,length(marketMinusRF));
dependentVar = reshape( excessReturn(13:end,:) , size(monthlyReturns(13:end,:),2),size(monthlyReturns(13:end,:),1)  );



%dependentVar = dependentVar(48:end,:);
%independentVar = independentVar(48:end,:);

betaCoeff = zeros( size(independentVar,2)-49,size(dependentVar,1  ));
betaStError = betaCoeff;

%treatNull = isnan(dependentVar) | isnan(independentVar);



    %for i=1:size(independentVar,1)
    for i=49:size(independentVar,2)
        for j=1:size(dependentVar,1)      
%       i=49;
%       j=1;
            %regResult = regstats( dependentVar( ~treatNull()  ,j), independentVar(~treatNull() ,j) );
            %regResult = regstats( dependentVar( i-48:i,j), independentVar(i-48:i,j) );
            regResult = regstats( dependentVar(j, i-48:i), independentVar(1, i-48:i) );
            
            betaCoeff(i-48,j) = regResult.tstat.beta(2);
            betaStError(i-48,j) = betaCoeff(i-48,j)/regResult.tstat.se(2);
        end
    end
    
save( 'Betas.mat','betaCoeff' );   
    
%end