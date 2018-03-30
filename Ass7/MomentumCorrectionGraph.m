%%
clear;
clc;
close all;


%This is the matlab file for the diagram
%Essentially, the new part starts from line 196


%Load monthly returns matrix
%load( 'monthlyReturns.mat' );
load('AllReturnsWithDate.mat');
load('SizeWithDates.mat');

    %Removing Last Month for 2018
    monthlyReturns = AllReturnsWithDate(1:size(AllReturnsWithDate,1)-1 ,:);
    %Correct data for NaN elements with value -99.99
    monthlyReturns( monthlyReturns==-99.99 ) =nan;
    monthlyReturns(:,2:end)=monthlyReturns(:,2:end)/100;
    %monthlyReturns=monthlyReturns(3:end,:);
    
    %Clear unneccesary variable
    clear monthlyReturns1;
    
%Removing Last Month for 2018   
sizeAssets = SizeWithDates(1:size(AllReturnsWithDate,1)-1,:);

 %Correct data for NaN elements with value -99.99
 sizeAssets ( sizeAssets==-99.99 ) =NaN;
 %Normalization
 for i=1:size(sizeAssets,1)
    sizeAssets(i,2:end) = sizeAssets(i,2:end)./max(sizeAssets(i,2:end));
 end
 

 
%Load FF factors
%The last month has been removed for 2018
load('marketMinusRF.mat');

load('SMB.mat');
    %High to low book to market
load('HML.mat');
load('RF.MAT');

%Define number of Assets and number of total months
numAssets = size(monthlyReturns,2)-1;
numMonths = size(monthlyReturns,1 );

%Rewritting Code for Momentum Technique
%Estimation of Compounded returns
%I treated null values as 0 in the sense that these assets, for this month,
%will neither be winners nor losers
compoundedReturns=ones( numMonths, numAssets+1) ;
compoundedReturns(:,1) = monthlyReturns(:,1);

for i=2:numAssets+1
    for j=2:numMonths
        compoundedReturns(j,i)=(1+monthlyReturns(j,i) )*compoundedReturns(j-1,i) ; 
    end
end

%plot(sizeAssets(3:end,2:end)./ sizeAssets(2:end -1,2:end), monthlyReturns(1:end - 2,2:end));
plot(compoundedReturns(:,2:end));


%%
%Creting Momentum Criterion Matrix
%At t I run Pt-1 / Pt-12
%In returns this means r_t-1 / r_t-11
momentumMatrixCriterion= nan( numMonths ,numAssets+1);
momentumMatrixCriterion(:,1)= monthlyReturns(:,1);
%Adding dates to Matrix
momentumMatrixCriterion( 12:end,2:end ) = -1 +compoundedReturns(12:end,2:end)./compoundedReturns(1:end-11,2:end);

%Skouras Code
% momentumMatrixCriterion=[nan(12,numAssets+1);...
%     -1+compoundedReturns(12:end-1,2:end)./compoundedReturns(1:end-12,2:end)];


%%
%Define number of winner and loser assets
numWinners=round(0.1*numAssets);
numLosers=round(0.1*numAssets);

%Correcting variables for dimention so as to be like momentumMatrix
marketMinusRF = marketMinusRF(1:numMonths);
SMB = SMB(1:numMonths);
HML = HML(1:numMonths);
RF = RF(1:numMonths);

%%

%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%

treatNanInSortCriterion = isnan(momentumMatrixCriterion);
treatNanInSortCriterionOrSize = isnan(momentumMatrixCriterion)| isnan(sizeAssets)  ;

%%%%%%%%% PART 3 %%%%%%%%%%%%%%%%%%%
%%%Panel A
%%I have already the returns of Momentum stored in matrix momentum

%Initialization of Return Matrices
ewR=nan(numMonths,2);
vwR=nan(numMonths,2);

%Adding dates to ewR and vwR matrices
ewR(:,1)= monthlyReturns(:,1);
vwR(:,1)= monthlyReturns(:,1);

%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 10;
numOfCategoriesCriterion2 = 5;

for i=12:numMonths-1
    %InvestmentDate=12;
    %Date for getting the return=13;
    %i=12;
    sortCriterion=momentumMatrixCriterion(i,2:end);
    % sortCriterion(isnan(sizeAssets(i,2:end)))=nan;
    sortCriterion(isnan(monthlyReturns(i+1,:)))=nan;

    [sortedMatrixofReturns,indexSorted] = sort(sortCriterion ,2);
      
    numOfAssetsOfMonth=sum(~isnan(sortedMatrixofReturns),2);
    numInDecileCriterion1=round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
    
    isLong=1 + indexSorted(...
          numOfAssetsOfMonth-numInDecileCriterion1+1:numOfAssetsOfMonth);
      
    isShort=1 + indexSorted(1:numInDecileCriterion1);
     
    ewR(i+1,2)...
        =nanmean(monthlyReturns(i+1,isLong))-...
         nanmean(monthlyReturns(i+1,isShort));  
     
     % equal weighted case
     
%      weightLong= ones(1,numOfAssetsOfMonth)/numOfAssetsOfMonth;
%      weightShort=ones(1,numOfAssetsOfMonth)/numOfAssetsOfMonth;
     
     
%       ewR(i+1,2)...
%         =sum(monthlyReturns(i+1,isLong).* weightLong)-...
%          sum(monthlyReturns(i+1,isShort).*weightShort); 
     

     weightLong= sizeAssets(i,isLong)/sum(sizeAssets(i,isLong));
     weightShort= sizeAssets(i,isShort)/sum(sizeAssets(i,isShort));
 
%      weightLong= sizeAssets(i+1,isLong)/sum(sizeAssets(i+1,isLong));
%      weightShort= sizeAssets(i+1,isShort)/sum(sizeAssets(i+1,isShort));
     
     
    vwR(i+1,2)...
        =sum(monthlyReturns(i+1,isLong).* weightLong)-...
         sum(monthlyReturns(i+1,isShort).*weightShort);      
                            
                     
end

ewRwithNoNan=ewR(:,2);ewRwithNoNan(isnan(ewRwithNoNan))=0;
vwRwithNoNan=vwR(:,2);vwRwithNoNan(isnan(vwRwithNoNan))=0;

figure(1)
plot(datenum(num2str(ewR(:,1)),'yyyymm'),cumprod(1+ewRwithNoNan));hold on
plot(datenum(num2str(vwR(:,1)),'yyyymm'),cumprod(1+vwRwithNoNan));
datetick
legend({ 'Equally Weighted Momentum Return','Value Weighted Momentum Return' } )
hold off


figure(2)
ewRExcessReturns = ewR(:,2) - (RF/100);ewRExcessReturns(isnan(ewRExcessReturns))=0;
vwRExcessReturns = vwR(:,2) - (RF/100);vwRExcessReturns(isnan(vwRExcessReturns))=0;


yyaxis left
plot(datenum(num2str(vwR(:,1)),'yyyymm'),cumprod(1+vwRExcessReturns) );hold on
ylabel('Cumulative MOM%')

yyaxis right
plot(datenum(num2str(vwR(:,1)),'yyyymm'),cumsum(log(1+vwRExcessReturns)));

xlabel('Date');
legend( {'Compounded excess return(left axis)' , 'Cumulative log excess return(right axis)'} );
datetick
hold off

%Save Equally weighted and Value weighted Momentum Returns
save( 'ewR.mat' , 'ewR' );
save( 'vwR.mat' , 'vwR' );




%Code for estimating Min, Max,Skewness & Kyrtosis
tempMomentumMatrixCriterion = momentumMatrixCriterion(:,2:end);
treatNan = isnan(tempMomentumMatrixCriterion);

%treatNan2 = isnan(momentumMatrix)| isnan(sizeAssets( 13:end ,:))  ;
averageValue = nanmean( nanmean(tempMomentumMatrixCriterion,2) );
minValue = nanmin(nanmin(tempMomentumMatrixCriterion));
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

for i = 12: numMonths
    %i=1;
      [sortedMatrixofReturns,indexMatrix] = sort( tempMomentumMatrixCriterion(i, ~treatNan(i,:)),2);
      numOfAssetsOfMonth = size(sortedMatrixofReturns,2);
     
      min5Value(1,i) = prctile(sortedMatrixofReturns,5);
      min25Value(1,i)= prctile(sortedMatrixofReturns,25);
      medianValue(1,i) = median( sortedMatrixofReturns);
      min75Value(1,i) = prctile(sortedMatrixofReturns,75);
      min95Value(1,i) =  prctile(sortedMatrixofReturns,95);
      maxValue(1,i) = max( sortedMatrixofReturns ); 
      numberOfAssets(1,i) = numOfAssetsOfMonth;
      stDeviation(1,i)= std( sortedMatrixofReturns );
      skewnessVector(1,i) = skewness(sortedMatrixofReturns);
      kurtosisVector(1,i) = kurtosis(sortedMatrixofReturns);
      
end

averageValue = averageValue*100
min5 = mean(min5Value)*100
min25 = mean(min25Value)*100
medianV = mean(medianValue)*100
min75 = mean(min75Value)*100
min95 = mean( min95Value)*100
maxV = mean(maxValue)*100
numberOfA = sum(numberOfAssets)
stDeviationA = mean(stDeviation)
skewnessA = mean(skewnessVector)
kurtosisA = mean(kurtosisVector)



%I have already the returns of Momentum stored in matrix momentum

[sortedMatrixofReturns,indexMatrix] = sort( tempMomentumMatrixCriterion,2);

%Initialization of Matrix for each percentile

mom =zeros(10,numMonths);

%Estimation of mean Realized Return for each Percentile
for i=1:numMonths
    
      [sortedMatrixofReturns,indexMatrix] = sort( tempMomentumMatrixCriterion(i, ~treatNan(i,:)),2);
      numOfAssetsOfMonth = size(sortedMatrixofReturns,2);
      mom(1,i)= nanmean( sortedMatrixofReturns(1: round(0.1*numOfAssetsOfMonth) ));
      mom(2,i)= nanmean( sortedMatrixofReturns(round(0.1*numOfAssetsOfMonth) +1 : round(0.2*numOfAssetsOfMonth) ));
      mom(3,i)= nanmean( sortedMatrixofReturns(round(0.2*numOfAssetsOfMonth)+1 : round(0.3*numOfAssetsOfMonth) ));
      mom(4,i)= nanmean( sortedMatrixofReturns(round(0.3*numOfAssetsOfMonth)+1: round(0.4*numOfAssetsOfMonth) ));
      mom(5,i)= nanmean( sortedMatrixofReturns(round(0.4*numOfAssetsOfMonth)+1: round(0.5*numOfAssetsOfMonth) ));
      mom(6,i)= nanmean( sortedMatrixofReturns(round(0.5*numOfAssetsOfMonth)+1: round(0.6*numOfAssetsOfMonth) ));
      mom(7,i)= nanmean( sortedMatrixofReturns(round(0.6*numOfAssetsOfMonth)+1: round(0.7*numOfAssetsOfMonth) ));
      mom(8,i)= nanmean( sortedMatrixofReturns(round(0.7*numOfAssetsOfMonth)+1: round(0.8*numOfAssetsOfMonth) ));
      mom(9,i)= nanmean( sortedMatrixofReturns(round(0.8*numOfAssetsOfMonth)+1: round(0.9*numOfAssetsOfMonth) ));
      mom(10,i)= nanmean( sortedMatrixofReturns(round(0.9*numOfAssetsOfMonth)+1: end ));
    
end

averageMOM = nanmean(mom,2)*100;
