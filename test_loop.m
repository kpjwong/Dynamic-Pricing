% test passenger inner loop
% say there are 2 options

% GMM
% load sim_data.mat;
% N = 100000; iter = 10;
% param_mat = zeros(7,iter);
% param = [10.750, -2.501, -11.908, -34.136, -2.873, -1.147, -7.783, -1.188, 0.549];
% 
% gen_data;
% 
% alpha = param(1:2); beta = param(3:4); gamma = param(5:6); delta = param(7:8); sigma = param(9);


for zz_idx = 14
    data = est_data; zz = z2zz(data(:,3),2);

    data = data(zz==zz_idx & data(:,2)==1,:);
    U = rand(size(data,1),2);
    z1 = data(:,3); z2 = data(:,4); taxi_P = data(:,7); taxi_T = data(:,8); uber_P = data(:,9); zz1 = z2zz(z1,2);
    dist = data(:,5); uber_T = data(:,10); mta_P = data(:,11); mta_T = data(:,12); zz2 = z2zz(z2,2);

    eps = zeros(1000,1);
    taxi_w = taxi_T-data(:,6); uber_w = uber_T-data(:,6);
    N = size(data,1);
    mta_w = zeros(N,1); mta_t = zeros(N,1);
    for n = 1:N
        mta_w(n) = wait_mta_zz(20*(data(n,2)-1) + zz1(n));
        mta_t(n) = mta_travel_time(z1(n),data(n,4));
    end

    eta = .01;
    
    for i = 1:999
        if i == 1
            D_taxi_tilde = N/3;
            D_uber_tilde = N/3;
        else
            D_taxi_tilde = eta*D_taxi + (1-eta)*D_taxi_tilde
            D_uber_tilde = eta*D_uber + (1-eta)*D_uber_tilde;
        end
        taxi_w_tilde = 10*60^(-0.3)*D_taxi_tilde^.4;
        uber_w_tilde = 10*40^(-0.8)*D_uber_tilde^.8;
        U_taxi = beta(1)*taxi_P + gamma(1)*data(:,6) + delta(1)*taxi_w_tilde;
        U_uber = beta(1)*uber_P + gamma(1)*data(:,6) + delta(1)*uber_w_tilde;
        U_nest = [U_taxi,U_uber]-max([U_taxi,U_uber],[],2);
        U_ride = alpha(1) + alpha(2) * (zz1 == 19 | zz2 == 19) + max([U_taxi,U_uber],[],2) + sigma(1)*log(sum(exp(U_nest/sigma(1)),2));
        ride_prob = exp(U_nest/sigma(1))./sum(exp(U_nest/sigma(1)),2);
        ride_choice = sum(bsxfun(@gt, U(:,1), cumsum(ride_prob,2)),2)+1;
        U_mta =  beta(2)*mta_P + gamma(2)*mta_t + delta(2)*mta_w;
        U_vect = [U_ride,U_mta]-max([U_ride,U_mta],[],2);
        nest_prob = exp(U_vect)./sum(exp(U_vect),2);
        nest_choice = sum(bsxfun(@gt, U(:,2), cumsum(nest_prob,2)),2)+1;
        sim_choice = zeros(N,1);
        sim_choice(nest_choice==2) = 3;
        sim_choice(nest_choice==1) = ride_choice(nest_choice==1);
        D_taxi = N*exp(U_taxi)/(exp(U_taxi)+exp(U_uber))*ride_prob;
        D_uber = N*;

        taxi_w = 10*60^(-0.3)*sum(sim_choice==1)^.4;
        uber_w = 10*40^(-0.8)*sum(sim_choice==2)^.8;

        eps(i+1) = norm([taxi_w_tilde-taxi_w,uber_w-uber_w_tilde]);

    end
    figure;
    title([num2str(zz_idx)])
    plot(eps)
end
% sim_mkt_share_taxi = zeros(19,24);
% sim_mkt_share_uber = zeros(19,24);
% for z_id = 1:19
%     for h = 0:23
%         sim_mkt_share_taxi(z_id,h+1) = sum(sim_choice==1 & zz==z_id & h_data==h)/size(data,1);
%         sim_mkr_share_uber(z_id,h+1) = sum(sim_choice==2 & zz==z_id & h_data==h)/size(data,1);
%     end
% end