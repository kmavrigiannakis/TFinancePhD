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





%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 10;
numOfCategoriesCriterion2 = 5;


%Initialization of Matrices
averageEwR = zeros(12,1);
averageVwR = averageEwR;
ValuesForTable = zeros(12,10);



for monthsAhead = 1:12
    %monthsAhead=1
    %Initialization of Return Matrices
    ewR=nan(numMonths,2);
    vwR=nan(numMonths,2);
    

    %Adding dates to ewR and vwR matrices
    ewR(:,1)= monthlyReturns(:,1);
    vwR(:,1)= monthlyReturns(:,1);
    
    
for i=12:numMonths - monthsAhead
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
        =nanmean(monthlyReturns(i+monthsAhead,isLong))-...
         nanmean(monthlyReturns(i+monthsAhead,isShort));  
     
     
%      ewR(i+1,2)...
%          =nanmean((compoundedReturns(i+monthsAhead,isLong)./compoundedReturns(i,isLong))...
%      -( compoundedReturns(i+ monthsAhead, isShort )./compoundedReturns(i,isShort) ));
     

     weightLong= sizeAssets(i,isLong)/sum(sizeAssets(i,isLong));
     weightShort= sizeAssets(i,isShort)/sum(sizeAssets(i,isShort));
 
%      weightLong= sizeAssets(i+1,isLong)/sum(sizeAssets(i+1,isLong));
%      weightShort= sizeAssets(i+1,isShort)/sum(sizeAssets(i+1,isShort));
     
     
    vwR(i+1,2)...
        =sum(monthlyReturns(i+monthsAhead,isLong).* weightLong)-...
         sum(monthlyReturns(i+monthsAhead,isShort).*weightShort);      

%     vwR(i+1,2)...
%         =sum(   (compoundedReturns(i+monthsAhead,isLong)./compoundedReturns(i,isLong) ).* weightLong)-...
%          sum( (compoundedReturns(i+ monthsAhead, isShort )./compoundedReturns(i,isShort) ).*weightShort);      

         
end

    averageEwR(monthsAhead) = nanmean(ewR(:,2))*100;
    averageVwR(monthsAhead) = nanmean(vwR(:,2))*100;
    
   
 %Regression Part

% ewRExcessReturns = ewR(:,2) - (RF/100);
% vwRExcessReturns = vwR(:,2) - (RF/100);

ewRExcessReturns = ewR(:,2)*100 - RF;
vwRExcessReturns = vwR(:,2)*100 - RF;
 
treatNaN = isnan(ewRExcessReturns);


 %Equally Weighted
 %CAPM
 [covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNaN),ewRExcessReturns(~treatNaN) );
 capmTstatisticsE = capmCoeffE(1)/capmNwErrorsE(1);
 
 %FAMA - MACBETH
 famaIndependentVar = [ marketMinusRF(~treatNaN), SMB(~treatNaN), HML(~treatNaN)  ];
 [covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar,ewRExcessReturns(~treatNaN) );
 famaTstatisticsE = famaCoeffE(1)/famaNwErrorsE(1);
 
 %Value Weighted
 %CAPM
 treatNaN = isnan(vwRExcessReturns);
 [covarCapm,capmNwErrorsV, capmCoeffV ] =  hac( marketMinusRF(~treatNaN),vwRExcessReturns(~treatNaN) );
 capmTstatisticsV = capmCoeffV(1)/capmNwErrorsV(1);
 
  %FAMA - MACBETH
 famaIndependentVar = [ marketMinusRF(~treatNaN), SMB(~treatNaN), HML(~treatNaN)  ];
 [covarFAMA,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar,vwRExcessReturns(~treatNaN) );
 famaTstatisticsV = famaCoeffV(1)/famaNwErrorsV(1);
 
       
 ValuesForTable(monthsAhead,1) = capmCoeffE(1);
 ValuesForTable(monthsAhead,2) = capmTstatisticsE;
 ValuesForTable(monthsAhead,3) = famaCoeffE(1);
 ValuesForTable(monthsAhead,4) = famaTstatisticsE;
 ValuesForTable(monthsAhead,5) = capmCoeffV(1);
 ValuesForTable(monthsAhead,6) = capmTstatisticsV;
 ValuesForTable(monthsAhead,7) = famaCoeffV(1);
 ValuesForTable(monthsAhead,8) = famaTstatisticsV;
 ValuesForTable(monthsAhead,9) = averageEwR(monthsAhead);
 ValuesForTable(monthsAhead,10) = averageVwR(monthsAhead);
    

end


% 
% figure(1)
% plot(datenum(num2str(ewR(:,1)),'yyyymm'),cumprod(1+ewRwithNoNan));hold on
% plot(datenum(num2str(vwR(:,1)),'yyyymm'),cumprod(1+vwRwithNoNan));
% datetick
% legend({ 'Equally Weighted Momentum Return','Value Weighted Momentum Return' } )
% hold off
% 


xlswrite('MonthsAhead',ValuesForTable);


