function CP = cal_AP(CLBP_SMCH,trainIDs, trainClassIDs, testIDs, testClassIDs)
% Authors: Zhenhua Guo, Lei Zhang, and David Zhang
    trains = CLBP_SMCH(trainIDs,:);
    tests = CLBP_SMCH(testIDs,:);
    trainNum = size(trains,1);
    testNum = size(tests,1);
    DM = zeros(testNum,trainNum);
    for i=1:testNum;
        test = tests(i,:);        
        DM(i,:) = distMATChiSquare(trains,test)';
    end
    CP=ClassifyOnNN(DM,trainClassIDs,testClassIDs);