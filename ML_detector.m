function [X_hat]= ML_detector(y,H,QAM_table, N)
y = y.';

n = N;            % Number of loops
v = ones (1, n);   % Index vector
Q = length(QAM_table);
Output = zeros (size(repmat (Q, n)));
ready = false;
metric = 100000;

while ~ ready
    index = sub2ind(size (Output), v);
    x_tmp = QAM_table(index);
    metric_tmp = sqrt((y - H*x_tmp.')'*(y - H*x_tmp.'));
    if metric_tmp<metric
        X_hat = x_tmp;
        metric = metric_tmp;
    end
    % Update the index vector:
    ready = true;       % Assume that the WHILE loop is ready
    for k = 1: n
        v (k) = v (k) + 1;
        if v (k) <= Q
            ready = false;
            break ;          % v (k) increased successfully, leave "for k" loop
        end
        v (k) = 1;         % v (k) reached the limit, reset it
    end
end

% metric = 100000;
% for l = 1:16
%     x_tmp(1) = QAM_table(l);
%     Esti_y1(:,1) = y - H(:,1).*x_tmp(1);
%     for m = 1:16
%         x_tmp(2) = QAM_table(m);
%         Esti_y2(:,1) = Esti_y1(:,1) - H(:,2).*x_tmp(2);
%         for n = 1:16
%             x_tmp(3) = QAM_table(n);
%             Esti_y3(:,1) = Esti_y2(:,1) - H(:,3).*x_tmp(3);
%             for o = 1:16
%                 x_tmp(4) = QAM_table(o);
%                 Esti_y4(:,1) = Esti_y3(:,1) - H(:,4).*x_tmp(4);
%                 metric_tmp = sqrt(Esti_y4(:,1)'*Esti_y4(:,1));
%                 metric_tmp2 = sqrt((y - H*x_tmp.')'*(y - H*x_tmp.'));
%                 if metric_tmp<metric
%                     X_hat = x_tmp; metric = metric_tmp;
%                 end
%             end
%         end
%     end
% end


