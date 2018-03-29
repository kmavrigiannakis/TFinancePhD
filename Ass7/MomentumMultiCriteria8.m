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


load('momentumMatrixCriterion.mat');


%Correcting variables for dimention so as to be like momentumMatrix
marketMinusRF = marketMinusRF(1:numMonths);
SMB = SMB(1:numMonths);
HML = HML(1:numMonths);
RF = RF(1:numMonths);


%Defining Sorting Criteria

%For Size
sortCriterion1 = sizeAssets(:,2:end);

%For BM
% load('BM.mat');
% BM = [ monthlyReturns(:,1) ,BM(1:1098,:)  ];
% sortCriterion1 = BM(:,2:end);

%For Betas
% load('Betas.mat');
% BetasNew =  nan( numMonths , numAssets+1);
% BetasNew(:,1) = monthlyReturns(:,1);
% BetasNew( size(monthlyReturns,1) -  size(betaCoeff,1) +1:end,2:end) = betaCoeff ;
% sortCriterion1 = BetasNew(:,2:end);

%Criterion 2 is always Momentum
sortCriterion2 = momentumMatrixCriterion(:,2:end);

%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 5;
numOfCategoriesCriterion2 = 5;

%%
%Main Part for the Estimation of Return
%For Depending Sorts


%Initialization of Return Matrices
ewR=nan(numMonths,2);
vwR=nan(numMonths,2);

%Adding dates to ewR and vwR matrices
ewR(:,1)= monthlyReturns(:,1);
vwR(:,1)= monthlyReturns(:,1);

%Initialization of Matrix of Returns
ewROfDecile = nan( numMonths, numOfCategoriesCriterion1+1,numOfCategoriesCriterion2+1);
vwROfDecile = ewROfDecile;

for i=12:numMonths-1
    %InvestmentDate=12;
    %Date for getting the return=13;
    %i=12;
  
    %Sorting according to Criterion 1
    [sortedMatrixofReturns,indexSortedCriterion1] = sort( sortCriterion1(i,:),2);
    
%     %Debug
%     sortedMatrixofReturns
%     indexSortedCriterion1
    
%     if (i==100)
%     a = sortedMatrixofReturns();
%     sum(and( ~isnan(sortCriterion1(i,:)) ,~isnan(sortCriterion2(i,:))))  
%     
%     end
    
    %Define how I will treat Nulls of two Criteria
    %numOfAssetsOfMonth = sum(~isnan(sortedMatrixofReturns),2);
    %numOfAssetsOfMonth = sum(and( ~isnan(sortCriterion1(i,:)) ,~isnan(sortCriterion2(i,:))))  ;
    numOfAssetsOfMonth = sum(or( ~isnan(sortCriterion1(i,:)), ~isnan(sortCriterion2(i,:))))  ;
 
%      numAssetsInDecileCriterion1 = round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
%      numAssetsInDecileCriterion2 = floor(numAssetsInDecileCriterion1 / numOfCategoriesCriterion2);
%     
%      %Debug
%      numOfAssetsOfMonth
%      numAssetsInDecileCriterion1
%      numAssetsInDecileCriterion2
     
     
%     numAssetsInDecileCriterion1 = floor(numOfAssetsOfMonth / numOfCategoriesCriterion1);
%     numAssetsInDecileCriterion2 = round(numAssetsInDecileCriterion1 / numOfCategoriesCriterion2);

    numAssetsInDecileCriterion1 = round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
    numAssetsInDecileCriterion2 = round(numAssetsInDecileCriterion1 / numOfCategoriesCriterion2);
%     
    
%     numAssetsInDecileCriterion1 = floor(numOfAssetsOfMonth / numOfCategoriesCriterion1);
%     numAssetsInDecileCriterion2 = floor(numAssetsInDecileCriterion1 / numOfCategoriesCriterion2);
%    
   
   
    %Sorting according to Criterion 2
    %Indexing for each Decile
    stepCriterion1 = numOfAssetsOfMonth;
    for j=numOfCategoriesCriterion1:-1:1
        
        %j=4;
%         if(j==4)
%             numOfAssetsOfMonth
%             numAssetsInDecileCriterion1
%             numAssetsInDecileCriterion2
%         end
        
        if (j~=1)
            [sortedMatrixofReturns,indexSortedCriterion2] = sort( sortCriterion2(i,indexSortedCriterion1( stepCriterion1 - numAssetsInDecileCriterion1+1: stepCriterion1 )  ),2);
        
%             %Debug
%             sortedMatrixofReturnsTemp(j,:)=sortedMatrixofReturns;
%             indexSortedCriterion2Temp(j,:)=indexSortedCriterion2;
%             
%             if i==100
%                 b = sortedMatrixofReturns
%                 numAssetsInDecileCriterion2
%             end
            stepCriterion1 = stepCriterion1 - numAssetsInDecileCriterion1;
        else
           
          [sortedMatrixofReturns,indexSortedCriterion2] = sort( sortCriterion2(i,indexSortedCriterion1( 1:numAssetsInDecileCriterion1  )  ),2);
        
%             %Debug
%             sortedMatrixofReturnsTemp(j,:)=sortedMatrixofReturns;
%             indexSortedCriterion2Temp(j,:)=indexSortedCriterion2;
%         
        
        end
    
    
    
   % stepCriterion2 = length(indexSortedCriterion2) - sum(isnan( sortedMatrixofReturns )  )   ;
        stepCriterion2 = length(indexSortedCriterion2)   ;  
        for k =numOfCategoriesCriterion2:-1:1
        
         if (k~=1)
             
%              stepCriterion2Initial = stepCriterion2
%              length(indexSortedCriterion2)
%              numAssetsInDecileCriterion2
             %ewROfDecile(i+1,k,j) = nanmean(monthlyReturns(i+1, 1+ indexSortedCriterion1( 1+ indexSortedCriterion2( stepCriterion2: stepCriterion2 + numAssetsInDecileCriterion2-1  )  ) ));
                ewROfDecile(i+1,k,j) = nanmean(monthlyReturns(i+1, ...
                    1+ indexSortedCriterion1( indexSortedCriterion2...
                    ( stepCriterion2 - numAssetsInDecileCriterion2+1: stepCriterion2  )  ) ));
                 
%                 stepCriterion2
%                 a=stepCriterion2 - numAssetsInDecileCriterion2+1
                
                weight = sizeAssets(i, 1+indexSortedCriterion1( indexSortedCriterion2...
                    ( stepCriterion2 - numAssetsInDecileCriterion2+1: stepCriterion2  ) ))...
                    /nansum(sizeAssets(i, 1+indexSortedCriterion1( indexSortedCriterion2...
                    ( stepCriterion2 - numAssetsInDecileCriterion2+1: stepCriterion2  )  ) )); 
                vwROfDecile(i+1,k,j) = nansum( monthlyReturns(i+1, 1+indexSortedCriterion1...
                    ( indexSortedCriterion2( stepCriterion2 - numAssetsInDecileCriterion2+1: stepCriterion2 )  ) )  .*weight );

                stepCriterion2 = stepCriterion2 - numAssetsInDecileCriterion2;
                
         else
             
                ewROfDecile(i+1,k,j) = nanmean(monthlyReturns(i+1,  1+indexSortedCriterion1( indexSortedCriterion2( 1:numAssetsInDecileCriterion2 ))));
               
                weight = sizeAssets(i, 1 + indexSortedCriterion1( indexSortedCriterion2...
                    ( 1:numAssetsInDecileCriterion2 ) ))...
                    /nansum(sizeAssets(i, 1+ indexSortedCriterion1( indexSortedCriterion2...
                    ( 1:numAssetsInDecileCriterion2  )))); 

                vwROfDecile(i+1,k,j) = nansum( monthlyReturns(i+1, 1 + indexSortedCriterion1...
                    ( indexSortedCriterion2( 1:numAssetsInDecileCriterion2 )  ) )  .*weight );
       
         end
            
             
        end
        
        
%         %Estimation of Momentum Portfolio
%         %1st way
          ewROfDecile(i+1,6,j) = ewROfDecile(i+1,5,j) - ewROfDecile(i+1,1,j);
          vwROfDecile(i+1,6,j) = vwROfDecile(i+1,5,j) - vwROfDecile(i+1,1,j);

          %2nd way
%         treatNan1 = sum(isnan( sortedMatrixofReturns));
%         isLong=1 + indexSortedCriterion1( indexSortedCriterion2...
%                (length(indexSortedCriterion2) - treatNan1 - numAssetsInDecileCriterion2+1: length(indexSortedCriterion2) - treatNan1 ));
%       
%         isShort=1 + indexSortedCriterion1( indexSortedCriterion2...
%                 (1:numAssetsInDecileCriterion2));
%      
% 
%         %Calculation for Equally weighted
%         ewROfDecile(i+1,6,j)...
%             =nanmean(monthlyReturns(i+1,isLong))-...
%              nanmean(monthlyReturns(i+1,isShort));  
%          
% %          a(i,:)=monthlyReturns(i+1,isLong);
% %          b(i,:)=monthlyReturns(i+1,isShort);
% %      
%         %Calculation for Value weighted
%         weightLong= sizeAssets(i,isLong)/sum(sizeAssets(i,isLong));
%         weightShort= sizeAssets(i,isShort)/sum(sizeAssets(i,isShort));
% 
%     
%         vwROfDecile(i+1,6,j)...
%         =sum(monthlyReturns(i+1,isLong).* weightLong)-...
%          sum(monthlyReturns(i+1,isShort).*weightShort);      
                            
    end    
    
%     ewROfDecile(i+1,:,6) = nanmean( ewROfDecile(i+1,:,:),2 );
%     squeeze( ewROfDecile(i+1,:,6) )
%     
    
    for m=1:6
        ewROfDecile(i+1,m,6) = nanmean( squeeze(ewROfDecile(i+1,m,1:5)) ,1 );
        vwROfDecile(i+1,m,6) = nanmean( squeeze(vwROfDecile(i+1,m,1:5)) ,1 );
    end
        
    
end
% 
% ewRwithNoNan=ewR(:,2);ewRwithNoNan(isnan(ewRwithNoNan))=0;
% vwRwithNoNan=vwR(:,2);vwRwithNoNan(isnan(vwRwithNoNan))=0;
% 
% figure(1)
% plot(datenum(num2str(ewR(:,1)),'yyyymm'),cumprod(1+ewRwithNoNan));hold on
% plot(datenum(num2str(vwR(:,1)),'yyyymm'),cumprod(1+vwRwithNoNan));
% datetick
% legend({ 'Equally Weighted Momentum Return','Value Weighted Momentum Return' } )
% hold off
% 
% 
%Ploting of Momentum Portfolio (for 5th decile)
ewRwithNoNan=ewROfDecile(:,6,5);ewRwithNoNan(isnan(ewRwithNoNan))=0;
plot( cumprod(1+ ewRwithNoNan ));

averageEwR = squeeze(nanmean( ewROfDecile,1 ));
%averageEwR(6,:) =  averageEwR(5,:) - averageEwR(1,:); 
%averageEwR(1:5,6) =  nanmean(averageEwR(1:5,1:5),2); 
%averageEwR(6,6) =  averageEwR(5,6) - averageEwR(1,6); 
averageEwR = averageEwR*100;


averageVwR = squeeze(nanmean( vwROfDecile,1 ));
%averageVwR(6,:) =  averageVwR(5,:) - averageVwR(1,:); 
%averageVwR(1:5,6) =  nanmean(averageVwR(1:5,1:5),2); 
%averageVwR(6,6) =  averageVwR(5,6) - averageVwR(1,6); 
averageVwR = averageVwR*100;


%% Write in Excel the Results
xlswrite('averageEwRSizeMomentum',averageEwR);
xlswrite('averageVwRSizeMomentum',averageVwR);

% xlswrite('averageEwR_BM_Momentum',averageEwR);
% xlswrite('averageVwR_BM_Momentum',averageVwR);

% xlswrite('averageEwR_Betas_Momentum',averageEwR);
% xlswrite('averageVwR_Betas_Momentum',averageVwR);
%%

% %%% Regression Part %%%%

%Initialization of Matrices
capmA = nan(1,numOfCategoriesCriterion1+1);
capmTstatisticsAlpha = capmA;
famaA =capmA;
famaTstatisticsAlpha = famaA;
capmTstatisticsAlphaV= famaA;
capmAV= famaA;
famaTstatisticsAlphaV= famaA;
famaAV= famaA;

% 
for j= 1:numOfCategoriesCriterion1+1
        
       %Treating Null in each iteration for appropriate regression
       treatNull = (isnan(ewROfDecile(:,6,j)));
       
       %Definition of independent variables for FAMA FRENCH Regression
       famaIndependentVar = [ marketMinusRF(~treatNull), SMB(~treatNull), HML(~treatNull)  ];
       
       %CAPM Regression for Equally Weighted Portfolio
       %[covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNull),  ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNull),  ewROfDecile(~treatNull,6,j)*100 );
       capmTstatisticsAlpha(j) = capmCoeffE(1)/capmNwErrorsE(1);
       capmA(j) = capmCoeffE(1);
       
       %Fama–MacBeth Regression for Equally Weighted Portfolio
       %[covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewROfDecile(~treatNull,6,j)*100   );
       famaTstatisticsAlpha(j) = famaCoeffE(1)/famaNwErrorsE(1);
       famaA(j) = famaCoeffE(1);
       
       %CAPM Regression for Value Weighted Portfolio
       %[covarCapmV,capmNwErrorsV, capmCoeffV ] =  hac( marketMinusRF(~treatNull),  vwROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [covarCapm,capmNwErrorsV, capmCoeffV ] =  hac( marketMinusRF(~treatNull),  vwROfDecile(~treatNull,6,j)*100 );
       capmTstatisticsAlphaV(j) = capmCoeffV(1)/capmNwErrorsV(1);
       capmAV(j) = capmCoeffV(1);
       
       %Fama–MacBeth Regression for  Value Weighted Portfolio
       %[covarFAMAV,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, vwROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [covarFAMA,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, vwROfDecile(~treatNull,6,j)*100   );
       famaTstatisticsAlphaV(j) = famaCoeffV(1)/famaNwErrorsV(1);
       famaAV(j) = famaCoeffV(1);
       
       
end

%regressResultsExcel =reshape( [ capmA ; capmTstatisticsAlpha; famaA; famaTstatisticsAlpha; capmAV;  capmTstatisticsAlphaV; famaAV; famaTstatisticsAlphaV ] ,6,8 );

%Write in Excel 
%%%For Size
%xlswrite('RegressionSizeAndMomentum',regressResultsExcel);

%For BM
%xlswrite('RegressionBMAndMomentum',regressResultsExcel);

%For Betas
%xlswrite('RegressionBetasAndMomentum',regressResultsExcel);