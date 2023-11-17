function plot_area_for_admissiable_poles(max_rise_time, max_overshoot, max_settlingtime)
% plots the area in which the poles are requiered to be in order to
% fullfil given requieremnts


w_n_min = 1.8/max_rise_time; % minimal distance from the origin
zeta_min = abs(log(max_overshoot))/sqrt((log(max_overshoot))^2 + pi^2); % minimal Damping
sigma_max = 4.6/max_settlingtime; % maximal real part

n_samples = 201;
figure();


% calculate bound according to maximal settling time
sigma_ts = ones(n_samples,1)'*(-sigma_max);
w_d_ts = linspace(-50,50,n_samples);

plot(sigma_ts,w_d_ts);
hold on

% calculate bound according to maximal overshoot
sigma_oS = linspace(-30,0,n_samples)*zeta_min;
w_d_oS = linspace(-30,0,n_samples)*sqrt(1-zeta_min.^2); 

plot([sigma_oS fliplr(sigma_oS)],[w_d_oS,-fliplr(w_d_oS)]);

% calculate bound according to maximal rise time
w_d_tr = linspace(-w_n_min,w_n_min,n_samples);
sigma_tr = -sqrt(w_n_min^2-w_d_tr.^2);

plot(sigma_tr,w_d_tr);


% Calculated the intersextion of the Graphs
sigma = sigma_oS(1);
w_d = w_d_oS(1);

[~,I] = min(abs(sigma_oS-sigma_ts));
[~,I2] = min(abs(sqrt(sigma_oS.^2+w_d_oS.^2)-sqrt(sigma_tr.^2+w_d_tr.^2)));
[min_distance3,I3] = min(abs(sigma_ts-sigma_tr));

% Make a case distingtion depenting on the position of the intersection
if sigma_oS(I) <= sigma_oS(I2)
    if(min_distance3 >= 0.01)
        sigma = [sigma sigma_oS(I)];
        w_d = [w_d w_d_oS(I)];
        sigma = [sigma sigma_oS(I)];
        w_d = [w_d -w_d_oS(I)];
    else
        sigma = [sigma sigma_oS(I) sigma_tr(I3)];
        w_d = [w_d w_d_oS(I) -w_d_tr(I3)];
        sigma = [sigma sigma_tr(length(sigma_tr)-(I3-1):(I3-1))];
        w_d = [w_d w_d_tr(length(w_d_tr)-(I3-1):(I3-1))]; 
        sigma = [sigma sigma_tr(I3) sigma_oS(I)];
        w_d = [w_d  w_d_tr(I3) -w_d_oS(I)];
    end
else
    sigma = [sigma sigma_oS(I2) ];
    w_d = [w_d w_d_oS(I2)];
    sigma = [sigma sigma_tr(length(sigma_tr)-(I2-1):(I2-1))];
    w_d = [w_d w_d_tr(length(w_d_tr)-(I2-1):(I2-1))]; 
    sigma = [sigma sigma_oS(I2) ];
    w_d = [w_d -w_d_oS(I2)];
end

sigma = [sigma sigma_oS(1)];
w_d = [w_d -w_d_oS(1)];

% Plot Area
fill(sigma,w_d,'b')
xlim([-10 0])
xlabel('Re')

ylim([-5 5])
ylabel('Im')

legend('settling time','maximal overshoot','rise time','Area of admissable poles')
end