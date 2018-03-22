clear;
clc;
close all;


%This is the matlab file for the diagram
%Essentially, the new part starts from line 196


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


%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%
%%%Panel A
%%I have already the returns of Momentum stored in matrix momentum

[sortedMatrixofReturns,indexMatrix] = sort( momentumMatrix,2);


%Estimation of realized Returns

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
      %[sortedMatrixofReturns,indexMatrixSize] = sort(montlyReturnsNoNullSize,2);
%       monthlyReturnsValueV = monthlyReturns( i+10, ~treatNan2(i,:) );
%       sizeAssetsV = sizeAssets( i+10, ~treatNan2(i,:) );
%       momentumMatrixV = momentumMatrix(i, ~treatNan2(i,:));
      
      
      monthlyReturnsValueV = monthlyReturns( i+10, ~treatNan2(i-1,:) );
      sizeAssetsV = sizeAssets( i+10, ~treatNan2(i-1,:) );
      momentumMatrixV = momentumMatrix(i-1, ~treatNan2(i-1,:));
      

      [sortedMatrixofReturns,indexMatrix] = sort( momentumMatrixV,2);
      
%       numOfNulls = 0;
%       temp= size(monthlyReturnsValueV,2);
%       while( isnan(sortedMatrixofReturns(temp))==1)
%           temp = temp -1 ;
%           numOfNulls = numOfNulls+1;
%       end
%       numOfAssetsOfMonthV =   numAssets - numOfNulls;
      
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
          realizedReturnValueW(temp,i-1) =  (sum( monthlyReturnsValueV(indexMatrix( newpositionV ) ) .* sizeAssetsV( indexMatrix(newpositionV ) ) ) )/ sum(sizeAssetsV( indexMatrix(newpositionV )))   ;
      end
      newpositionE= round(0.9 *numOfAssetsOfMonth) +1 :  numOfAssetsOfMonth ;
      realizedReturnEqualW(10,i-1) =  sum( monthlyReturns( i+10, indexMatrixE( newpositionE ) ) ) ;
      realizedReturnValueW(10,i-1) =  (sum( monthlyReturnsValueV(indexMatrix( newpositionV ) ) .* sizeAssetsV( indexMatrix(newpositionV ) ) ) )/ sum(sizeAssetsV( indexMatrix(newpositionV )))   ;

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
 
%%Essentially this is the new Part
%Plot Cumulative Excess Returns and Log Excess Returns
 
 %cumReturns = cumprod((realizedReturnValueW(11,:)/10 )+1)-1 ;
 figure(1)
 
 excessReturns = realizedReturnEqualW(11,:) - transpose(RF);
 cumReturns = cumprod( (excessReturns/100) +1  )-1 ;
 cumLogReturns = cumsum( log( 1+ excessReturns) ) ;
 
 cumTime =1: size(cumReturns,2);
 
 %Setting Time
 startingDate = datetime('30-June-1927');
 endingDate = datetime(2017,12,31);
 matrixOfDates = linspace( startingDate, endingDate, size(cumReturns,2)  );
 
 hold on
 
 yyaxis left
 plot(matrixOfDates,cumReturns);  
 ylabel('Cumulative MOM%')
 
 
 % plot(cumTime,cumLogReturns);  
 yyaxis right
 plot(matrixOfDates,real(cumLogReturns));  
 legend( {'Compounded excess return' , 'Cumulative log excess return'   } );
 
 xlabel('Date')
 
hold off
 
 
