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
load('BM.mat');
BM = [ monthlyReturns(:,1) ,BM(1:1098,:)  ];
%sortCriterion1 = BM(:,2:end);

%For Betas
  load('Betas.mat');
  BetasNew =  nan( numMonths , numAssets+1);
  BetasNew(:,1) = monthlyReturns(:,1);
  BetasNew( size(monthlyReturns,1) -  size(betaCoeff,1) +1:end,2:end) = betaCoeff ;
 %sortCriterion1 = BetasNew(:,2:end);

 
%%%%%%%%%% CORRELATIONS %%%%%%%%%%%%%%%%%%
%Timeseries Averages of Cross Sectional Correlation
corrMomentumSize = nanmean(nanmean(corr(momentumMatrixCriterion( 12:end,2:end )  ,sizeAssets(12:end,2:end) )  )   );

% b =  nanmean(nanmean(corr( reshape(  momentumMatrixCriterion( 12:end,2:end ), size(momentumMatrixCriterion,2)-1 , ...
%     size(momentumMatrixCriterion,1)-11 ),...
%     reshape( sizeAssets(12:end,2:end), size(momentumMatrixCriterion,2)-1 , size(momentumMatrixCriterion,1)-11  )) ));

corrMomentumSizeSpearman = nanmean(nanmean(corr(momentumMatrixCriterion( 12:end,2:end ), sizeAssets(12:end,2:end)  ,'type','Spearman') ));


corrMomentumBetas = nanmean(nanmean(corr(momentumMatrixCriterion( 61:end,2:end ), BetasNew(61:end,2:end) )  )   );
corrMomentumBetasSpearman = nanmean(nanmean(corr(momentumMatrixCriterion( 61:end,2:end ), BetasNew(61:end,2:end)  ,'type','Spearman') ));

corrMomentumBM = nanmean(nanmean(corr(momentumMatrixCriterion( 12:end,2:end ), BM(12:end,2:end) )  )   );
corrMomentumBMSpearman = nanmean(nanmean(corr(momentumMatrixCriterion( 12:end,2:end ), BM(12:end,2:end)  ,'type','Spearman') ));

% 
% corrMomentumSize = corr(nanmean(momentumMatrixCriterion( 12:end,2:end ),2),nanmean(sizeAssets(12:end,2:end),2));
% corrMomentumSizeSpearman = corr(nanmean(momentumMatrixCriterion( 12:end,2:end ),2), nanmean(sizeAssets(12:end,2:end),2) ,'type','Spearman'); 
% 
% corrMomentumBetas = corr(nanmean(momentumMatrixCriterion( 61:end,2:end ),2),nanmean( BetasNew(61:end,2:end),2));
% corrMomentumBetasSpearman = corr(nanmean(momentumMatrixCriterion( 61:end,2:end ),2), nanmean( BetasNew(61:end,2:end),2) ,'type','Spearman'); 
% 
% corrMomentumBM = corr(nanmean(momentumMatrixCriterion( 12:end,2:end ),2),nanmean(BM(12:end,2:end),2));
% corrMomentumBMSpearman = corr(nanmean(momentumMatrixCriterion( 12:end,2:end ),2), nanmean(BM(12:end,2:end),2)  ,'type','Spearman'); 
%  
 
 
%Criterion 2 is always Momentum
sortCriterion2 = momentumMatrixCriterion(:,2:end);

%Define Number of Categories for each Criterion
numOfCategoriesCriterion1 = 3;
numOfCategoriesCriterion2 = 3;


%Main Part for the Estimation of Return
%For InDependent Sorts


%Initialization of Return Matrices
ewR=nan(numMonths,4);
vwR=nan(numMonths,4);


%Initialization of Matrix of Returns of 9*9 matrix
ewROfQuintile = nan( numMonths, numOfCategoriesCriterion1+1);
vwROfQuintile = ewROfQuintile;

for i=12:numMonths-1
    %InvestmentDate=12;
    %Date for getting the return=13;
    %i=12;
  
    %Define how I will treat Nulls of two Criteria
    %numOfAssetsOfMonth = sum(~isnan(sortedMatrixofReturns),2);
    %numOfAssetsOfMonth = sum(and( ~isnan(sortCriterion1(i,:)) ,~isnan(sortCriterion2(i,:))))  ;
    %numOfAssetsOfMonth = sum(or( ~isnan(sortCriterion1(i,:)), ~isnan(sortCriterion2(i,:))))  ;
    numOfAssetsOfMonth = sum( ~isnan(sortCriterion1(i,:)) & ...
        ~isnan(sortCriterion2(i,:)) & ~isnan(monthlyReturns( i+1, 2:end) )) ;  
     
   numAssetsInQuantileCriterion1 = round(numOfAssetsOfMonth / numOfCategoriesCriterion1);
   numAssetsInQuantileCriterion2 = round(numOfAssetsOfMonth / numOfCategoriesCriterion2);

   
    %Sorting according to Criterion 1
    [sortedMatrixofReturnsCriterion1,indexSortedCriterion1] = sort( sortCriterion1(i,:),2);
    [sortedMatrixofReturnsCriterion2,indexSortedCriterion2] = sort( sortCriterion2(i,:),2);
       
    %Sorting according to Criterion 2
    %Indexing for each Decile
    stepCriterion1 = numOfAssetsOfMonth;
    
    for j1 = numOfCategoriesCriterion1:-1:1
            
        if (j1~=1)
            
            positionInIndex1 = stepCriterion1 - numAssetsInQuantileCriterion1+1: stepCriterion1;
            
        else

            positionInIndex1 =  1:numAssetsInQuantileCriterion1;
            
        end
        stepCriterion2 = numOfAssetsOfMonth;
        for j2 = numOfCategoriesCriterion2:-1:1
            %Estimation of Returns for Equally Weighted Portfolio
%             ewROfDecile(i+1,j) = nanmean(monthlyReturns(i+1, ...
%                     1+  indexAssetsOfCategory( positionInIndex )));
%              
            if(j2~=1)
                
            positionInIndex2 = stepCriterion2 - numAssetsInQuantileCriterion2+1: stepCriterion2;
                
            else
                
            positionInIndex2 =  1:numAssetsInQuantileCriterion2;
            
            end
            
            indexAssetsOfCategory = intersect( indexSortedCriterion1( positionInIndex1 ) , indexSortedCriterion2(positionInIndex2) );
            
            ewROfQuintile(j1,j2) = nanmean(monthlyReturns(i+1, ...
                    1+  indexAssetsOfCategory));
            
            %Estimation of Returns for Value Weighted Portfolio
%             weight = sizeAssets(i, 1 + indexAssetsOfCategory( positionInIndex ))...
%                     /nansum(sizeAssets(i, 1+ indexAssetsOfCategory( positionInIndex ))) ; 
% 
%             vwROfDecile(i+1,k,j) = nansum( monthlyReturns(i+1, 1 +  indexAssetsOfCategory( positionInIndex )).*weight );
%                 
            
            weight = sizeAssets(i, 1 + indexAssetsOfCategory)...
                    /nansum(sizeAssets(i, 1+ indexAssetsOfCategory )) ; 

            vwROfQuintile(j1,j2) = nansum( monthlyReturns(i+1, 1 + ...
                indexAssetsOfCategory).*weight );
                
            stepCriterion2 = stepCriterion2 - numAssetsInQuantileCriterion2;
        
%         ewROfQuintile(i+1,numOfCategoriesCriterion1+1) = nanmean( ewROfQuintile(i+1,1:numOfCategoriesCriterion1) ,2 );
%         vwROfQuintile(i+1,numOfCategoriesCriterion1+1) = nanmean( vwROfQuintile(i+1,1:numOfCategoriesCriterion1) ,2 );
%     
        end
        
        stepCriterion1 = stepCriterion1 - numAssetsInQuantileCriterion1;
    
         ewR(i,j1) = ewROfQuintile(j1,3) - ewROfQuintile(j1,1);
         vwR(i,j1) = vwROfQuintile(j1,3) - vwROfQuintile(j1,1);

        
    end
        
    ewR(i,4) = nanmean(ewR(i,1:3));
    vwR(i,4) = nanmean(vwR(i,1:3));
   
    
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
% %Ploting of Momentum Portfolio (for 5th decile)
% ewRwithNoNan=ewROfDecile(:,numOfCategoriesCriterion2+1,numOfCategoriesCriterion1);ewRwithNoNan(isnan(ewRwithNoNan))=0;
% plot( cumprod(1+ ewRwithNoNan ));
% 
averageEwR = squeeze(nanmean( ewR,1 ));
%averageEwR(6,:) =  averageEwR(5,:) - averageEwR(1,:); 
%averageEwR(1:5,6) =  nanmean(averageEwR(1:5,1:5),2); 
%averageEwR(6,6) =  averageEwR(5,6) - averageEwR(1,6); 
averageEwR = averageEwR*100;


averageVwR = squeeze(nanmean( vwR,1 ));
%averageVwR(6,:) =  averageVwR(5,:) - averageVwR(1,:); 
%averageVwR(1:5,6) =  nanmean(averageVwR(1:5,1:5),2); 
%averageVwR(6,6) =  averageVwR(5,6) - averageVwR(1,6); 
averageVwR = averageVwR*100;
% 
% 
% %% Write in Excel the Results
% xlswrite('averageEwRSizeMomentum',averageEwR);
% xlswrite('averageVwRSizeMomentum',averageVwR);
% 
% % xlswrite('averageEwR_BM_Momentum',averageEwR);
% % xlswrite('averageVwR_BM_Momentum',averageVwR);
% 
% % xlswrite('averageEwR_Betas_Momentum',averageEwR);
% % xlswrite('averageVwR_Betas_Momentum',averageVwR);
% %%
% 
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
onesTstatistics = nan(2,numOfCategoriesCriterion1+1);
onesRegression = ones( numMonths,1 );

% 
for j= 1:numOfCategoriesCriterion1+1
        
       %Treating Null in each iteration for appropriate regression
       treatNull = (isnan(ewR(:,j)));
       
       %Definition of independent variables for FAMA FRENCH Regression
       famaIndependentVar = [ marketMinusRF(~treatNull), SMB(~treatNull), HML(~treatNull)  ];
       
       
       %Regression with one's  
       [~,onesNwErrorsE, onesCoeffE ] =  hac( onesRegression(~treatNull),  ewR(~treatNull,j)*100, 'intercept', false );
       onesTstatistics(1,j) = onesCoeffE(1)/onesNwErrorsE(1);
       [~,onesNwErrorsV, onesCoeffV ] =  hac( onesRegression(~treatNull),  vwR(~treatNull,j)*100, 'intercept', false );
       onesTstatistics(2,j) = onesCoeffV(1)/onesNwErrorsV(1);
       
       
       %CAPM Regression for Equally Weighted Portfolio
       %[covarCapm,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNull),  ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [~,capmNwErrorsE, capmCoeffE ] =  hac( marketMinusRF(~treatNull),  ewR(~treatNull,j)*100 );
       capmTstatisticsAlpha(j) = capmCoeffE(1)/capmNwErrorsE(1);
       capmA(j) = capmCoeffE(1);
       
       %Fama–MacBeth Regression for Equally Weighted Portfolio
       %[covarFAMA,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [~,famaNwErrorsE, famaCoeffE ] =  hac( famaIndependentVar, ewR(~treatNull,j)*100   );
       famaTstatisticsAlpha(j) = famaCoeffE(1)/famaNwErrorsE(1);
       famaA(j) = famaCoeffE(1);
       
       %CAPM Regression for Value Weighted Portfolio
       %[covarCapmV,capmNwErrorsV, capmCoeffV ] =  hac( marketMinusRF(~treatNull),  vwROfDecile(~treatNull,6,j)*100 - RF(~treatNull) );
       [~,capmNwErrorsV, capmCoeffV ] =  hac( marketMinusRF(~treatNull),  vwR(~treatNull,j)*100 );
       capmTstatisticsAlphaV(j) = capmCoeffV(1)/capmNwErrorsV(1);
       capmAV(j) = capmCoeffV(1);
       
       %Fama–MacBeth Regression for  Value Weighted Portfolio
       %[covarFAMAV,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, vwROfDecile(~treatNull,6,j)*100 - RF(~treatNull)  );
       [~,famaNwErrorsV, famaCoeffV ] =  hac( famaIndependentVar, vwR(~treatNull,j)*100   );
       famaTstatisticsAlphaV(j) = famaCoeffV(1)/famaNwErrorsV(1);
       famaAV(j) = famaCoeffV(1);
       
       
end


%Write Data to Excel
forExcel = [averageEwR ; onesTstatistics(1,:) ; capmA ; capmTstatisticsAlpha ; famaA ; famaTstatisticsAlpha;  averageVwR ; onesTstatistics(2,:) ; capmAV ; capmTstatisticsAlphaV ; famaAV ; famaTstatisticsAlphaV;   ];
xlswrite('forExcel',forExcel);


%FAMA - MacBeth Regressions

% numOfIndependent = 2;
% NwErrors = nan( numMonths-11,numOfIndependent );
% Coeff = NwErrors;
% 
% for year = 12: numMonths
%     
%     IndependentVar = momentumMatrixCriterion(year,2:end);
%     
%     [~,NwErrors(year-11,:), Coeff(year-11,:) ] =  hac( IndependentVar,monthlyReturns(year,2:end) );
%     
% end    
% 
% 
% onesMatrix = ones(numMonths-11);
% for i=1:numOfIndependent
%     [~,onesNwErrors, onesCoeff ] =  hac( onesMatrix,  Coeff(:,i), 'intercept', false );
%     regressionTstatistics = onesCoeff./onesNwErrors;
% end
% 
% 

