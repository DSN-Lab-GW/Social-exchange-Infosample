function NLL = mymodelUncertainty_ORA(pars, data)
% beta1      =  pars(1);
% k          =  pars(2);
% alphaPrior =  pars(3);
% betaPrior  =  pars(4);
beta1      =  exp(pars(1));
k          =  1/(1+exp(-pars(2))); 
alphaPrior =  1+exp(pars(3)); 
betaPrior  =  1+exp(pars(4));
% alphaPrior =  1; 
% betaPrior  =  1;

% NLL = 0;
m   = 2;
T   = 25;

        DeltaQ = computeUncertainty_ORA(T, m, k, alphaPrior, betaPrior);
        
        trialidx   = find(data.red + data.green < T);
        
        thistime   = data.red(trialidx) + data.green(trialidx) + 1;
        thischoice = data.choice(trialidx);
        
        linearidx           = sub2ind(size(DeltaQ), data.green(trialidx) + 1, thistime);
        DeltaQ_vectorized   = DeltaQ(:);
        
        % Log likelihood
        prediction = 1./(1+exp(- thischoice .* (beta1 * DeltaQ_vectorized(linearidx))));
%         NLL = NLL - sum(log(prediction));
        NLL = sum(log(prediction+eps));

        
%     lowLimits =  [0   0 1  1 ];
%     highLimits = [inf 1 inf inf];
