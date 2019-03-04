% clear;
% load sim_data.mat
uber_tr_tinvariant(59,:) = uber_tr_tinvariant(58,:);
% generate dmd data for taxi and uber, day 1, zz = [1:5,7:16]
% data format [choice t z1 z2 dist dur t_price t_time u_price u_time m_price m_time]
zz_taxi = z2zz(taxi_dmd_data(:,3),2);
est_taxi_dmd_data = taxi_dmd_data(taxi_dmd_data(:,1)==1 & zz_taxi~=6 & zz_taxi < 17, 2:end);
est_taxi_dmd_data_boro = taxi_dmd_data(taxi_dmd_data(:,1)==1 & (zz_taxi==6 | zz_taxi >= 17) & zz_taxi < 20, 2:end);
est_taxi_dmd_data = [ones(size(est_taxi_dmd_data,1),1), est_taxi_dmd_data, zeros(size(est_taxi_dmd_data,1),5)];
est_taxi_dmd_data_boro = [ones(size(est_taxi_dmd_data_boro,1),1), est_taxi_dmd_data_boro, zeros(size(est_taxi_dmd_data_boro,1),5)];

% data append before sampling
zz_uber = z2zz(uber_dmd_data(:,3),2);
est_uber_dmd_data = uber_dmd_data(uber_dmd_data(:,1)==1 & zz_uber~=6 & zz_uber < 17, 2:end);
est_uber_dmd_data_boro = uber_dmd_data(uber_dmd_data(:,1)==1 & (zz_uber==6 | zz_uber >= 17) & zz_uber < 20, 2:end);
est_uber_dmd_data = [2*ones(size(est_uber_dmd_data,1),1), est_uber_dmd_data, zeros(size(est_uber_dmd_data,1),9)];
est_uber_dmd_data_boro = [2*ones(size(est_uber_dmd_data_boro,1),1), est_uber_dmd_data_boro, zeros(size(est_uber_dmd_data_boro,1),9)];

zz_mta = z2zz(mta_est_data(:,2),2);
est_mta_dmd_data = mta_est_data(zz_mta~=6 & zz_mta < 17,:);
est_mta_dmd_data_boro = mta_est_data((zz_mta==6 | zz_mta >= 17) & zz_mta < 20,:);
est_mta_dmd_data = [3*ones(size(est_mta_dmd_data,1),1), est_mta_dmd_data, zeros(size(est_mta_dmd_data,1),9)];
est_mta_dmd_data_boro = [3*ones(size(est_mta_dmd_data_boro,1),1), est_mta_dmd_data_boro, zeros(size(est_mta_dmd_data_boro,1),9)]; 

% sample N passengers
est_data = [est_taxi_dmd_data;est_uber_dmd_data;est_mta_dmd_data];
est_data_boro = [est_taxi_dmd_data_boro;est_uber_dmd_data_boro;est_mta_dmd_data_boro];

est_data_all = [[est_data,ones(size(est_data,1),1)]; [est_data_boro,zeros(size(est_data_boro,1),1)]];
est_data_all = est_data_all(est_data_all(:,3)>1 & z_b(est_data_all(:,3))<5,:);
M = size(est_data_all,1); id = randi([1,M],N,1); est_data_all = est_data_all(id,:);
est_data = est_data_all(est_data_all(:,end)==1,1:end-1);
est_data_boro = est_data_all(est_data_all(:,end)==0,1:end-1);

t = est_data(:,2); h = floor((t-1)/6) + 1; z1 = est_data(:,3);
zz = z2zz(est_data(:,3),2);
for i = 1:size(est_data,1)
    if est_data(i,1)~=1
        if est_data(i,1)==2
            % draw destinations for uber trips
            est_data(i,4) = sum(rand() > cumsum(uber_tr_h(z1(i),:,h(i)))) + 1;
            if est_data(i,4) > 263
                est_data(i,4) = sum(rand() > cumsum(uber_tr_tinvariant(z1(i),:))) + 1;
            end
        elseif est_data(i,1)==3
            % draw destinations for mta trips
            est_data(i,4) = sum(rand() > cumsum(mta_tr(z1(i),:,h(i)))) + 1;
        end
        z2 = est_data(i,4);
        % draw trip distances and trip durations for uber trips
        if mean_trip_distance(z1(i),z2,h(i)) > 0
            est_data(i,5) = max(normrnd(mean_trip_distance(z1(i),z2,h(i)), sd_trip_distance(z1(i),z2,h(i))),0.3);
        else
            est_data(i,5) = max(normrnd(distance_matrix(z1(i),z2),1),0.3);
        end
        est_data(i,5) = min(est_data(i,5),60);
        if mean_trip_duration(z1(i),z2,h(i)) > 0
            est_data(i,6) = max(normrnd(mean_trip_duration(z1(i),z2,h(i)), sd_trip_duration(z1(i),z2,h(i))),3);
        else
            est_data(i,6) = max(60*est_data(i,5)/v_h(zz(i),h(i)),3);
        end
        % impute taxi price and travel time for taxi trips
        if mean_taxi_fare(z1(i),z2,h(i)) > 0
            est_data(i,7) = mean_taxi_fare(z1(i),z2,h(i));
        else
            est_data(i,7) = 1.18*(2.5 + 2.5 * est_data(i,5));
        end
    end
    z2 = est_data(i,4);
    % impute taxi travel time
    est_data(i,8) = exprnd(mean_taxi_wait(zz(i),h(i))) + est_data(i,6);
    % impute uber price and travel time for taxi trips
    est_data(i,9) = max(sm_uber(1,t(i),z1(i)) * (2.55 + 0.3 * est_data(i,6) + 1.75 * est_data(i,5)), 8);
    est_data(i,10) = wait_uber(1,t(i),z1(i)) + est_data(i,6)/pool_ratio(z1(i),h(i));
    % impute mta price and travel time
    est_data(i,11) = 2.75;
    est_data(i,12) = mta_walking_time(z1(i)) + wait_mta_zz(20*(t(i)-1) + zz(i)) + mta_travel_time(z1(i),z2);
end

% data format [choice t z1 z2 dist dur t_price t_time u_price u_time m_price m_time]
t = est_data_boro(:,2); h = floor((t-1)/6) + 1; z1 = est_data_boro(:,3); 
z1 = est_data_boro(:,3);
b = z_b(z1);
zz = z2zz(z1,2);
for i = 1:size(est_data_boro,1)
    % draw destinations for uber and mta trips
    if est_data_boro(i,1)~=1
        if est_data_boro(i,1)==2
            % draw destinations for uber trips
            est_data_boro(i,4) = sum(rand() > cumsum(uber_tr_h(z1(i),:,h(i)))) + 1;
            if est_data_boro(i,4) > 263
                est_data_boro(i,4) = sum(rand() > cumsum(uber_tr_tinvariant(z1(i),:))) + 1;
            end
        elseif est_data_boro(i,1)==3
            % draw destinations for mta trips
            est_data_boro(i,4) = sum(rand() > cumsum(mta_tr(z1(i),:,h(i)))) + 1;
        end
        z2 = est_data_boro(i,4);
        % draw trip distances and trip durations for uber trips
        if mean_trip_distance(z1(i),z2,h(i)) > 0
            est_data_boro(i,5) = max(normrnd(mean_trip_distance(z1(i),z2,h(i)), sd_trip_distance(z1(i),z2,h(i))),0.3);
        else
            est_data_boro(i,5) = max(normrnd(distance_matrix(z1(i),z2),1),0.3);
        end
        est_data_boro(i,5) = min(est_data_boro(i,5),60);
        if mean_trip_duration(z1(i),z2,h(i)) > 0
            est_data_boro(i,6) = max(normrnd(mean_trip_duration(z1(i),z2,h(i)), sd_trip_duration(z1(i),z2,h(i))),3);
        else
            if zz(i) < 19
                v = v_h(zz(i),h(i));
            else
                v = vb_h(b(i),h(i));
            end
            est_data_boro(i,6) = max(60*est_data_boro(i,5)/v,3);
        end
        % impute taxi price and travel time for taxi trips
        if mean_taxi_fare(z1(i),z2,h(i)) > 0
            est_data_boro(i,7) = mean_taxi_fare(z1(i),z2,h(i));
        else
            est_data_boro(i,7) = 1.18*(2.5 + 2.5 * est_data_boro(i,5));
        end
    end
    if w_z(z1(i),h(i)) > 0
        est_data_boro(i,8) = exprnd(w_z(z1(i),h(i))) + est_data_boro(i,6);
    else
        if zz(i) < 19
            est_data_boro(i,8) = exprnd(w_zz(zz(i),h(i))) + est_data_boro(i,6);
        else
            est_data_boro(i,8) = exprnd(w_b(b(i),h(i))) + est_data_boro(i,6);
        end
    end
    est_data_boro(i,9) = max(sm_uber(1,t(i),z1(i)) * (2.55 + 0.3 * est_data_boro(i,6) + 1.75 * est_data_boro(i,5)), 8);
    est_data_boro(i,10) = wait_uber(1,t(i),z1(i)) + est_data_boro(i,6)/pool_ratio(z1(i),h(i));
    % impute mta price and travel time
    est_data_boro(i,11) = 2.75;
    est_data_boro(i,12) = mta_walking_time(z1(i)) + wait_mta_zz(20*(t(i)-1) + zz(i)) + mta_travel_time(z1(i),z2);
end

est_data = [est_data;est_data_boro];