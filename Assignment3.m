load('size.mat');
sizeIndustry= Size;
clear Size;

load('booktomarket.mat');

load('DATAFinanceExercise2.mat');
monthlyReturns = AverageValueWeightedReturnsMonthly;
clear AverageValueWeightedReturnsMonthly;


numWinners=round(0.1*size(monthlyReturns,2));
numLosers=round(0.1*size(monthlyReturns,2));
%Equally Weighted Portfolio
weight = 1/ numWinners;

%Estimation of STReversal
% sortedMatrixofReturns = zeros(size(monthlyReturns));
% indexMatrix=zeros(size(monthlyReturns));

% for time=2:size(monthlyReturns,1)
%     [sortedMatrixofReturns(time,:),indexMatrix(time,:)] =sort( monthlyReturns(time,:));
% end

[sortedMatrixofReturns,indexMatrix] = sort( monthlyReturns,2);
numAssets = size(monthlyReturns,2);


PorReturn= zeros(1, size(monthlyReturns,1 )-1 );
for time =2: size(monthlyReturns,1)
  %time=2;
% %time=1;
%     %Remove of Nan Elements for each month
      tempNan = isnan(sortedMatrixofReturns(time-1,:)); 
%     %tempVector = sortedMatrixofReturns(time, 1:( numAssets - sum(tempNan,2) ));
%     
%         %Estimation of Return for Losers
%         PorReturn(time)= sum(weight*tempVector(1:numLosers));
%         %Estimation of Return for Winers
%         PorReturn(time)= PorReturn -sum(weight*tempVector(1, length(tempVector)-numWinners+1 : length(tempVector)  ));
%     %end

        tempPorReturn=0;
        for i=1:numLosers
            tempPorReturn = tempPorReturn  + weight* monthlyReturns(time,indexMatrix(time-1,i));
        end
        
        for i=  numAssets - (sum(tempNan,2)) - numWinners +1 :  numAssets - (sum(tempNan,2))
            tempPorReturn = tempPorReturn-weight* monthlyReturns(time,indexMatrix(time-1,i));
        end

        PorReturn(time-1)= tempPorReturn;
end

figure(1)
correctTime =2: size(monthlyReturns,1);
plot( log(correctTime),cumsum(PorReturn));  
%cdfplot(PorReturn )
datetick
axis tight

%plot(cumsum(PorReturn));  


datetick

    %
%  for numAsset =1: (numWinners+numLosers)
%    
%     while(numAsset <= numLosers )
%             
%     end
%     
% end
% 



