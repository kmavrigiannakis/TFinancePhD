clear;
clc;
close all;

%% Initialization %%%
%%Loading of appropriate Variables
%Load monthly returns matrix
load( 'monthlyReturns.mat' );
    %Remove data for 2018
    monthlyReturns = monthlyReturns1(1:size(monthlyReturns1,1)-1 ,:);
    %Correct data for NaN elements with value -99.99
    monthlyReturns( monthlyReturns==-99.99 ) =0;
    monthlyReturns=monthlyReturns/100;
    %Clear unneccesary variable
    clear monthlyReturns1;
load('size.mat');  
sizeAssets = Size(1:1098,:);
 %Correct data for NaN elements with value -99.99
 sizeAssets ( sizeAssets==-99.99 ) =NaN;
 for i=1:size(sizeAssets,1)
    sizeAssets(i,:) = sizeAssets(i,:)./max(sizeAssets(i,:));
 end
 
clear Size;

%Loading returns from the Momentum results of the previous assignment
load('momentumMatrix.mat');
%Load FF factors
load('marketMinusRF.mat');

load('SMB.mat');
    %High to low book to market
load('HML.mat');
load('RF.MAT');

%Esimating the cumulative returns of momentum technique if needed
cumMomentum = cumprod(momentum+1)-1;

%Define number of Assets and number of total months
numAssets = size(monthlyReturns,2);
numMonths = size(monthlyReturns,1 );

%Rewritting Code for Momentum Technique
%Estimation of Compounded returns
%I treated null values as 0 in the sense that these assets, for this month,
%will neither be winners nor losers
compoundedReturns=ones( size(monthlyReturns,1 ), size(monthlyReturns,2)) ;
for i=1:size(monthlyReturns,2 )
    for j=2:size(monthlyReturns,1 )
%         if j==2
%             prices(j-1,i)
%         end
        compoundedReturns(j,i)=(1+monthlyReturns(j,i) )*compoundedReturns(j,i) ; 
    end
end

%Estimation of Pt-1 / Pt-12
momentumMatrix=zeros(size(monthlyReturns,1 )-12 , size(monthlyReturns,2 ));
for currentAssetNum = 1:numAssets
    for i = 13: numMonths
        momentumMatrix(i-12,currentAssetNum) =  compoundedReturns(i-1,currentAssetNum)/ compoundedReturns(i-12,currentAssetNum);
    end
end
%Correcting matrix so as to remove later the nan elements
%momentumMatrix( momentumMatrix ==1 ) = -100;
momentumMatrix( momentumMatrix ==1 ) = NaN;


%Define number of winner and loser assets
numWinners=round(0.1*numAssets);
numLosers=round(0.1*numAssets);

%Correcting variables for dimention so as to be like momentumMatrix
marketMinusRF = marketMinusRF(1:size(momentumMatrix,1));
SMB = SMB(1:size(momentumMatrix,1));
HML = HML(1:size(momentumMatrix,1));
RF = RF(1:size(momentumMatrix,1));


%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%
%Code for estimating Min, Max,Skewness & Kyrtosis
treatNan = isnan(momentumMatrix);
treatNan2 = isnan(momentumMatrix)| isnan(sizeAssets( 13:end ,:))  ;
averageValue = mean( nanmean(momentumMatrix,2) );
minValue = min(min(momentumMatrix));
%[sortedMatrixofReturns,indexMatrix] = sort( momentumMatrix,2);

%Initialization of Necessary Vectors
min5Value = zeros(1,numMonths);
min25Value = min5Value ;
medianValue = min5Value ;
min75Value = min5Value ;
min95Value = min5Value ;
maxValue = min5Value ;
numberOfAssets = min5Value ;
stDeviation =  min5Value ;
skewnessVector =  min5Value ;
kurtosisVector = min5Value ;

%Estimation of min values of Formating Period

for i = 1: size(momentumMatrix,1)
    %i=1;
      [sortedMatrixofReturns,indexMatrix] = sort( momentumMatrix(i, ~treatNan(i,:)),2);
      numOfAssetsOfMonth = size(sortedMatrixofReturns,2);
     
      min5Value(1,i) = prctile(sortedMatrixofReturns,5);
      min25Value(1,i)= prctile(sortedMatrixofReturns,25);
      medianValue(1,i) = median( sortedMatrixofReturns(1: round(0.25*numOfAssetsOfMonth) ));
      min75Value(1,i) = prctile(sortedMatrixofReturns,75);
      min95Value(1,i) =  prctile(sortedMatrixofReturns,95);
      maxValue(1,i) = max( sortedMatrixofReturns ); 
      numberOfAssets(1,i) = numOfAssetsOfMonth;
      stDeviation(1,i)= std( sortedMatrixofReturns );
      skewnessVector(1,i) = skewness(sortedMatrixofReturns);
      kurtosisVector(1,i) = kurtosis(sortedMatrixofReturns);
      
end

min5 = mean(min5Value)
min25 = mean(min25Value)
medianV = mean(medianValue) 
min75 = mean(min75Value)
min95 = mean( min95Value)
maxV = mean(maxValue)
numberOfA = sum(numberOfAssets)
stDeviationA = mean(stDeviation)
skewnessA = mean(skewnessVector)
kurtosisA = mean(kurtosisVector)

%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%
corrMomentumSize = corr(nanmean(momentumMatrix,2), SMB );
corrMomentumSizeSpearman = corr(nanmean(momentumMatrix,2), SMB  ,'type','Spearman');



%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%
%%%Panel A
%%I have already the returns of Momentum stored in matrix momentum

%Ploting the cumulative returns
% figure(1)
% correctTime =2: size(momentumMatrix,1);
% plot( (correctTime),cumprod(momentum+1)-1 ,'r'); 

[sortedMatrixofReturns,indexMatrix] = sort( momentumMatrix,2);

%Initialization of Matrix for each percentile

mom =zeros(10,size(momentumMatrix,1));

%Estimation of mean Realized Return for each Percentile
for i=1:size(momentumMatrix,1)
    
      [sortedMatrixofReturns,indexMatrix] = sort( momentumMatrix(i, ~treatNan(i,:)),2);
      numOfAssetsOfMonth = size(sortedMatrixofReturns,2);
      mom(1,i)= mean( sortedMatrixofReturns(1: round(0.1*numOfAssetsOfMonth) ));
      mom(2,i)= mean( sortedMatrixofReturns(round(0.1*numOfAssetsOfMonth) +1 : round(0.2*numOfAssetsOfMonth) ));
      mom(3,i)= mean( sortedMatrixofReturns(round(0.2*numOfAssetsOfMonth)+1 : round(0.3*numOfAssetsOfMonth) ));
      mom(4,i)= mean( sortedMatrixofReturns(round(0.3*numOfAssetsOfMonth)+1: round(0.4*numOfAssetsOfMonth) ));
      mom(5,i)= mean( sortedMatrixofReturns(round(0.4*numOfAssetsOfMonth)+1: round(0.5*numOfAssetsOfMonth) ));
      mom(6,i)= mean( sortedMatrixofReturns(round(0.5*numOfAssetsOfMonth)+1: round(0.6*numOfAssetsOfMonth) ));
      mom(7,i)= mean( sortedMatrixofReturns(round(0.6*numOfAssetsOfMonth)+1: round(0.7*numOfAssetsOfMonth) ));
      mom(8,i)= mean( sortedMatrixofReturns(round(0.7*numOfAssetsOfMonth)+1: round(0.8*numOfAssetsOfMonth) ));
      mom(9,i)= mean( sortedMatrixofReturns(round(0.8*numOfAssetsOfMonth)+1: round(0.9*numOfAssetsOfMonth) ));
      mom(10,i)= mean( sortedMatrixofReturns(round(0.9*numOfAssetsOfMonth)+1: end ));
    
end

averageMOM = mean(mom,2);

%Panel B

realizedReturnEqualW = zeros(11,size(momentumMatrix,1));
realizedReturnValueW = realizedReturnEqualW;

%New Returns for removing Null elements from size

for i=2:size(momentumMatrix,1)
    %i=2;
      [sortedMatrixofReturns,indexMatrixE] = sort( momentumMatrix(i-1,:),2);
      
      %Treating for Null for Equally Weighted
      numOfNulls = 0;
      temp= numAssets;
      while( isnan(sortedMatrixofReturns(temp))==1)
          temp = temp -1 ;
          numOfNulls = numOfNulls+1;
      end
      numOfAssetsOfMonth =   numAssets - numOfNulls;  
      
      %Treating for Null for Value Weighted
      %Code for synchronizing momentum matrix, size matrix, and monthly returns
      monthlyReturnsValueV = monthlyReturns( i+10, ~treatNan2(i-1,:) );
      sizeAssetsV = sizeAssets( i+10, ~treatNan2(i-1,:) );
      momentumMatrixV = momentumMatrix(i-1, ~treatNan2(i-1,:));
      [sortedMatrixofReturns,indexMatrix] = sort( momentumMatrixV,2);
      numOfAssetsOfMonthV = size(monthlyReturnsValueV,2);
      %The starting position for matrix monthly returns is at 12
      %Because the criterion is Pt-1 / Pt-12
      %I begin from i=2 so I need to add 10
      %Equally Weighted Portfolio
      
      temp=0;
      for j = 0:0.1:0.8
          temp = temp+1;
          newpositionE= round(j*numOfAssetsOfMonth) +1 : round( (j+0.1)*numOfAssetsOfMonth);
          newpositionV= round(j*numOfAssetsOfMonthV) +1 : round( (j+0.1)*numOfAssetsOfMonthV);
          realizedReturnEqualW(temp,i-1) =  sum( monthlyReturns( i+10, indexMatrixE(newpositionE)));
          realizedReturnValueW(temp,i-1) =  10*(sum( monthlyReturnsValueV(indexMatrix( newpositionV ) ) .* sizeAssetsV( indexMatrix(newpositionV ) ) ) )/ sum(sizeAssetsV( indexMatrix(newpositionV )))   ;
      end
      newpositionE= round(0.9 *numOfAssetsOfMonth) +1 :  numOfAssetsOfMonth ;
      realizedReturnEqualW(10,i-1) =  sum( monthlyReturns( i+10, indexMatrixE( newpositionE ) ) ) ;
      realizedReturnValueW(10,i-1) =  10*(sum( monthlyReturnsValueV(indexMatrix( newpositionV ) ) .* sizeAssetsV( indexMatrix(newpositionV ) ) ) )/ sum(sizeAssetsV( indexMatrix(newpositionV )))   ;

      %Essentially my return using the Momentum Technique
      realizedReturnValueW(11,i-1) =  realizedReturnValueW(10,i-1) - realizedReturnValueW(1,i-1);
      realizedReturnEqualW(11,i-1) =  realizedReturnEqualW(10,i-1) - realizedReturnEqualW(1,i-1);
     
end

averageRealizedReturnE = mean(realizedReturnEqualW,2);
averageRealizedReturnW = nanmean(realizedReturnValueW,2);


%Regression Part

 %Equally Weighted
 %CAPM
 capmRegE = regstats( realizedReturnEqualW(11,:) , marketMinusRF );
 capmCoeffE = capmRegE.tstat.beta;
 residualsE = capmRegE.r;
 capmNwErrorsE = nwse(residualsE,marketMinusRF );
 capmTstatisticsE = capmCoeffE./capmNwErrorsE;
 
 %FAMA - MACBETH
 famaIndependentVar = [ marketMinusRF, SMB, HML  ];
 famaRegE = regstats( realizedReturnEqualW(11,:) , famaIndependentVar );
 famaCoeffE = famaRegE.tstat.beta;
 residualsE = famaRegE.r;
 famaNwErrorsE = nwse(residualsE,famaIndependentVar );
 famaTstatisticsE = famaCoeffE./famaNwErrorsE;
 
 %Value Weighted
 %CAPM
 capmRegV = regstats( realizedReturnValueW(11,:) , marketMinusRF );
 capmCoeffV = capmRegV.tstat.beta;
 residualsV = capmRegV.r;
 capmNwErrorsV = nwse(residualsV,marketMinusRF );
 capmTstatisticsV = capmCoeffV./capmNwErrorsV;
 
 %FAMA - MACBETH
 famaRegV = regstats( realizedReturnValueW(11,:) , famaIndependentVar );
 famaCoeffV = famaRegV.tstat.beta;
 residualsV = famaRegV.r;
 famaNwErrorsV = nwse(residualsV,famaIndependentVar );
 famaTstatisticsV = famaCoeffV./famaNwErrorsV;
 
 
