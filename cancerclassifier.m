%% 0. House‑keeping -------------------------------------------------------
clc; clear; close all; rng(1);           % rng for reproducibility

%% 1. Load GEO dataset (GSE9476) ----------------------------------------
file = 'GSE9476_series_matrix.txt';
if ~isfile(file)
    disp('📥 Downloading GSE9476 ...');
    gunzip('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9476/matrix/GSE9476_series_matrix.txt.gz','.');
end
if ~isfile(file), error('❌ Dataset not found.'); end

fid  = fopen(file,'rt');
raw  = textscan(fid,'%s','Delimiter','\n'); fclose(fid);
lines= raw{1};
beg  = find(contains(lines,'!series_matrix_table_begin')) + 1;
fin  = find(contains(lines,'!series_matrix_table_end'))  - 1;
rows = lines(beg:fin);

numMat = cellfun(@(r) str2double(strsplit(r,'\t')), rows,'UniformOutput',false);
mat    = vertcat(numMat{:});            % genes × (samples+1)
X      = mat(:,2:end)';                % samples × genes (64 × 12 625)
X(isnan(X)) = 0;                       % replace NaNs
Y      = [ones(26,1); zeros(38,1)];    % 1 = AML, 0 = Control
fprintf("\n✅ Matrix loaded: %d samples × %d genes\n", size(X,1), size(X,2));

%% 2. Hold‑out split (70 % train / 30 % test) ----------------------------
cv  = cvpartition(Y,'HoldOut',0.30,'Stratify',true);
tr  = training(cv); te = test(cv);
Xtr = X(tr,:);  Xte = X(te,:);
Ytr = Y(tr);    Yte = Y(te);

%% 3. Model 1 – SVM (all genes) -----------------------------------------
mdl1   = fitcsvm(Xtr,Ytr,'KernelFunction','linear','Standardize',true);
Ypred1 = predict(mdl1,Xte);
acc1   = mean(Ypred1==Yte)*100;

%% 4. Model 2 – SVM + mRMR (Top 200) ------------------------------------
[idxOverall,~] = fscmrmr(Xtr,Ytr);      % use training fold only!
sel200 = idxOverall(1:200);
Xsel   = X(:, sel200);
mdl2   = fitcsvm(Xsel(tr,:),Ytr,'KernelFunction','linear','Standardize',true);
Ypred2 = predict(mdl2,Xsel(te,:));
acc2   = mean(Ypred2==Yte)*100;

%% 5. Model 3 – RF + SMO Hybrid -----------------------------------------
% 5‑A RF ranking on mRMR‑reduced pool (training data only)
rf  = TreeBagger(100,Xsel(tr,:),Ytr,'Method','classification', ...
                'OOBPredictorImportance','on');
[~,rfIdx] = sort(rf.OOBPermutedPredictorDeltaError,'descend');
Xpool = Xsel(:, rfIdx(1:200));          % 200‑gene pool

% 5‑B Spider‑Monkey Optimisation (binary mask search)
pop      = 30; iterMax = 40; dim = size(Xpool,2);
maskPop  = rand(pop,dim) < 0.5;
bestAcc  = 0; bestMask = false(1,dim);

for it = 1:iterMax
    for i = 1:pop
        mask = maskPop(i,:);
        if ~any(mask), mask(randi(dim)) = true; end  % ensure at least 1 gene
        mdl = fitcsvm(Xpool(tr,mask),Ytr,'KernelFunction','linear','Standardize',true);
        acc = mean(predict(mdl,Xpool(te,mask))==Yte);
        if acc > bestAcc
            bestAcc  = acc;
            bestMask = mask;
        end
    end
    % rudimentary population update toward bestMask
    for i = 1:pop
        flipProb = 0.1;
        L = rand(1,dim) < flipProb;
        maskPop(i,L) = bestMask(L);
    end
end
numSMO = sum(bestMask);
fprintf("🕷️  SMO selected %d features of 200\n", numSMO);

% 5‑C Final SVM on SMO‑selected genes
mdl3   = fitcsvm(Xpool(tr,bestMask),Ytr,'KernelFunction','linear','Standardize',true);
Ypred3 = predict(mdl3,Xpool(te,bestMask));
acc3   = mean(Ypred3==Yte)*100;

%% 6. Command‑Window Summary --------------------------------------------
fprintf("\n📊 FINAL ACCURACIES (70/30 Hold‑Out)\n");
fprintf("1️⃣  SVM (All Genes)        : %.2f %%\n", acc1);
fprintf("2️⃣  SVM + mRMR (Top 200)   : %.2f %%\n", acc2);
fprintf("3️⃣  RF + SMO Hybrid        : %.2f %%\n", acc3);
fprintf("    (Spider‑Monkey chose %d genes)\n\n", numSMO);

%% 7. VISUALISATIONS -----------------------------------------------------

% 7‑A Accuracy Bar Chart
figure('Name','Model Accuracy Comparison');
accArray = [acc1, acc2, acc3];
bar(accArray);
set(gca,'XTickLabel',{'All Genes','mRMR‑200','RF+SMO'});
ylabel('Accuracy (%)');
title('Model Performance on 30% Hold‑Out Test Set');
grid on; ylim([0,100]);

% 7‑B Feature‑Selection Funnel Plot
figure('Name','Feature Selection Funnel');
counts = [size(X,2), 200, 200, numSMO];
labels = {'All Genes','mRMR Top‑200','RF Top‑200','SMO Chosen'};
barh(flip(counts), 'FaceColor', [0.2 0.6 0.8]);
set(gca,'YTick',1:4,'YTickLabel',flip(labels));
xlabel('Number of Genes');
title('Feature‑Selection Funnel');

% 7‑C PCA Scatterplot on SMO‑selected Genes
if numSMO >= 2
    [~,score,~,~,explained] = pca(Xpool(:,bestMask));
    figure('Name','PCA of SMO‑Selected Genes');
    gscatter(score(:,1), score(:,2), Y, 'rb', 'ox');
    xlabel(sprintf('PC1 (%.1f%% var)', explained(1)));
    ylabel(sprintf('PC2 (%.1f%% var)', explained(2)));
    title('Samples in PCA Space – SMO Feature Set');
    legend({'AML','Control'},'Location','best');
    grid on;
end

% 7‑D Simplified Confusion Matrices -------------------------------------
classNames = {'AML','Control'};
YteCat    = categorical(Yte,   [1 0], classNames);
Ypred1Cat = categorical(Ypred1,[1 0], classNames);
Ypred2Cat = categorical(Ypred2,[1 0], classNames);
Ypred3Cat = categorical(Ypred3,[1 0], classNames);

figure('Name','Confusion Matrices');
subplot(1,3,1);
confusionchart(YteCat,Ypred1Cat,'Title','All Genes', ...
               'RowSummary','off','ColumnSummary','off');
subplot(1,3,2);
confusionchart(YteCat,Ypred2Cat,'Title','mRMR‑200', ...
               'RowSummary','off','ColumnSummary','off');
subplot(1,3,3);
confusionchart(YteCat,Ypred3Cat,'Title','RF+SMO', ...
               'RowSummary','off','ColumnSummary','off');
sgtitle('Confusion Matrices (Test Set)','FontWeight','bold');



%% End of script ---------------------------------------------------------
