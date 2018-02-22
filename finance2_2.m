clear;
clc;
close all;

load( 'finance2data.mat');
data= AverageValueWeightedReturnsMonthly;

load('IndustriesNames.mat');

clear AverageValueWeightedReturnsMonthly ;

Average = nanmean(data);
VarCovar = cov(data,'partialrows');


% Estimating Efficient Frontier with Default Constraints
p = Portfolio('AssetMean', Average, 'AssetCovar', VarCovar, 'LowerBound', 0, 'UpperBound',1, 'Budget', 1);

% Estimating Efficient Frontier with No Constraints
%p = Portfolio('AssetMean', Average, 'AssetCovar', VarCovar, 'LowerBound', -1, 'UpperBound',2, 'Budget', 1);

NumPor=100;
EffFrontier = estimateFrontier(p,NumPor);

figure(1);
plotFrontier(p, NumPor);
[risks , returns] = plotFrontier(p, NumPor);

figure(2);
plotFrontier(p, NumPor);
%Estimating the optimum Sharpe Ratio which is tangled with efficient
%frontier line
weights = estimateMaxSharpeRatio(p);
[risk, ret] = estimatePortMoments(p, weights);
hold on
plot(risk,ret,'*r');


%Estimating weights for each industry for each month which max Sharpe Ratio
MatrixMeanReturn = zeros( 1,size(data,1));

%
for i=1:size(data,1)
   %tempMin = nanmean( data(i,:))
   tempMin =2;
   for j=1:size(data,2)
      if ( (data(i,j)>0 )&& (data(i,j)<tempMin) )
         tempMin = data(i,j); 
      end
   MatrixMeanReturn(1,i)= tempMin;    
   end 
end

RiskFree = mean(MatrixMeanReturn);
SharpeRatioWeights = sharpe(data,RiskFree);

figure(3)
%pie(SharpeRatioWeights,pieOnesZeros,cellstr(IndustriesNames));
ShareRweightsPer=strings(1,49);
for i=1:1:49
ShareRweightsPer(1,i) = num2str(round(SharpeRatioWeights(1,i)*100,2));
end

combinedtxt = strcat( cellstr(IndustriesNames),ShareRweightsPer );

%pie(SharpeRatioWeights  ,pieOnesZeros, cellstr(combinedtxt));


pie(SharpeRatioWeights  ,cellstr(combinedtxt));


TableOfWeights= [ reshape( IndustriesNames, 49,1), round( reshape( SharpeRatioWeights, 49,1) ,2 )];

filename = 'C:\Users\kostmavr1986\Google Drive\PHD\Finance\Assignments\2\Weights.xlsx';
xlswrite(filename,TableOfWeights);

