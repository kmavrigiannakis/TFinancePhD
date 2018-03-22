%function [ averageRealizedReturnE, realizedReturnEqualW] = sortCriterion2E( momentumMatrix, monthlyReturns, criterion1, criterion2, numPercentiles1, numPercentiles2  )


%%%%%
clear;
clc;
close all;

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

load('betas.mat');

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


%Definining Criteria
criterion1 = sizeAssets;
%criterion1 = betaCoeff;

criterion2 = momentumMatrix;
numPercentiles1=5;
numPercentiles2=5;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%







load('momentumMatrix.mat');

treatNan = isnan(momentumMatrix);
%treatNan2 = isnan(momentumMatrix)| isnan(sizeAssets( 13:end ,:))  ;


%Initialize realized returns Matrix according to argument numPercentiles1
%Which is the number of categories I want
realizedReturnEqualW = zeros(size(momentumMatrix,1), numPercentiles1+1, numPercentiles2+1 );


%New Returns for removing Null elements from size
numAssets = size(momentumMatrix,2);

for i=2: min ( size(momentumMatrix,1), size(criterion1,1) )
   % i=2;
      [sortedMatrixofReturns,indexMatrixECriterion1] = sort( criterion1(i-1,:),2);
      
      %Treating for Null for Equally Weighted
      numOfNulls = sum(treatNan(i,:));
      numOfAssetsOfMonth =   numAssets - numOfNulls;  
           
      %The starting position for matrix monthly returns is at 12
      %Because the criterion is Pt-1 / Pt-12
      %I begin from i=2 so I need to add 10
      %Equally Weighted Portfolio
      
      temp1=0;
      numberOfAssetsinCategory1 = floor(numOfAssetsOfMonth/numPercentiles1);
      numberOfAssetsinCategory2 = floor(numberOfAssetsinCategory1/numPercentiles2);
      temp2 = 0;
      %for j = 1:numberOfAssetsinCategory1: (numOfAssetsOfMonth - numberOfAssetsinCategory1)
      %for j=1:numPercentiles2-1
      for j=1:numPercentiles2
        %j=  numPercentiles2;
          %temp1 = temp1+1;
          newpositionE= temp1+1: temp1 + numberOfAssetsinCategory1;
          temp1=temp1 + numberOfAssetsinCategory1;
          %newpositionE
          %newpositionE= j: numberOfAssetsinCategory1+j;
         % newpositionE
             
           [sortedCriterion2E,indexMatrixE] = sort(criterion2(i-1, newpositionE),2);
          % [sortedCriterion2E,indexMatrixE] = sort(criterion2(i-1, indexMatrixECriterion1),2);
         
          %for x=0:numberOfAssetsinCategory2:newpositionE
           for x=1:numPercentiles2-1
   
             %temp2=temp2 + numberOfAssetsinCategory2;
%            criterion2( i-1, newpositionE(indexMatrixE) )
%            monthlyReturns( i+10,
          
                %newpositionECriterion2 = x+1 : numberOfAssetsinCategory2+x;
                newpositionECriterion2 = temp2+1: temp2 + numberOfAssetsinCategory2;
                temp2=temp2 + numberOfAssetsinCategory2;
               % newpositionECriterion2
                
                %Estimation of Realized return for each subcategory
                %realizedReturnEqualW(i-1, j+1,x+1) =  sum( monthlyReturns( i+10, indexMatrixE(newpositionE)));
                %realizedReturnEqualW(i-1,temp1,temp2) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
                
                %realizedReturnEqualW(i-1,j,x) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
                realizedReturnEqualW(i-1,x,j) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
          
           end
           %temp2
           %temp2+1
           %length(newpositionE)
           %temp2+1
           %length(newpositionE)
           newpositionECriterion2 = temp2+1 : length(newpositionE) ;
          % newpositionECriterion2
           % temp2=temp2 + length(newpositionE);
           % temp2 = newpositionE(length(newpositionE));
          
            temp2 = 0;
            
           
            %realizedReturnEqualW(i-1, j,numPercentiles2 ) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
            realizedReturnEqualW(i-1, numPercentiles2,j ) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
           
  
      end
      %realizedReturnEqualW(i-1, numPercentiles1, numPercentiles2 ) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
      realizedReturnEqualW(i-1, numPercentiles2, numPercentiles1 ) =  sum( monthlyReturns( i+10,( indexMatrixECriterion1( newpositionE( indexMatrixE(newpositionECriterion2 ) ) )  ) ) );
       
      % realizedReturnEqualW(i-1, numPercentiles1, numPercentiles2 )
%       if( isnan(sum(monthlyReturns( i+10, indexMatrix( newposition ))) )==1 )
%          i 
%       end
      
%       if( isnan(sum(sizeAssets( i+10, indexMatrix(newpositionV ))))  ==1 )
%          i
%       end
%       
  
      %Calculation of last Column
      realizedReturnEqualW(i-1,1:numPercentiles2,numPercentiles1+1) = mean(realizedReturnEqualW(i-1, 1:numPercentiles2,1:numPercentiles1) ,3 );
      
      %Essentially my return using the Momentum Technique (10-1)
      realizedReturnEqualW(i-1, numPercentiles2+1,:) =  realizedReturnEqualW(i-1,numPercentiles2,:) - realizedReturnEqualW(i-1,1,:);
     
end

averageRealizedReturnE = squeeze(mean(realizedReturnEqualW,1));
 

%%% Regression Part

capmA = zeros(numPercentiles2+1,numPercentiles1+1);
capmTstatisticsAlpha = capmA;
famaA =capmA;
famaTstatisticsAlpha = famaA;

% if size(realizedReturnEqualW,1) < size(marketMinusRF,1)
%         
%         marketMinusRF = marketMinusRF(    )
% end

famaIndependentVar = [ marketMinusRF, SMB, HML  ];

for i=1: numPercentiles2+1
    for j= 1:numPercentiles1+1
        
%        [capmA(i,j),capmTstatisticsAlpha(i,j)] = capmRegstats(realizedReturnEqualW(11:end,i,j));
%        [famaA(i,j),famaTstatisticsAlpha(i,j)] = famaRegstats(realizedReturnEqualW(11:end,i,j));
%        
       [covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF,  realizedReturnEqualW(:,i,j) );
       capmTstatisticsAlpha(i,j) = capmCoeffE(1)/capmNwErrorsE(1);
       capmA(i,j) = capmCoeffE(1);
       
       [covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, realizedReturnEqualW(:,i,j));
       famaTstatisticsAlpha(i,j) = famaCoeffE(1)/famaNwErrorsE(1);
       famaA(i,j) = famaCoeffE(1);
       
    end
end


%Write in Excel
     forExcel = [ capmA(6,:) ; capmTstatisticsAlpha(6,:); famaA(6,:) ; famaTstatisticsAlpha(6,:)   ];
xlswrite('Table11_4_5',forExcel );

%Write tables to Excel
xlswrite('AverageReturns',averageRealizedReturnE);





