the_path =strcat(data_path,'\','*.csv'); %The path to the .csv files is passed to the MATLAB script via the labview program

cd(data_path);
dd = dir('*.csv');                           
                     

fileNames = {dd.name}; 

for i=1:length(fileNames)  %Changing comma to dot because German computer!

	file    = memmapfile( char(fileNames(i)), 'writable', true );
	comma   = uint8(',');
	point   = uint8('.');
	file.Data( transpose( file.Data==comma) ) = point;
	
	
	tab   = uint8('	');
	semicolon   = uint8(';');
	file.Data( transpose( file.Data==tab) ) = semicolon;
	
	
	new_line   = 10;
	space   = uint8(' ');
	file.Data( transpose( file.Data==new_line) ) = space;
	

end



data = cell(numel(fileNames),2);
data(:,1) = regexprep(fileNames, '.csv','');

for ii = 1:numel(fileNames)    
   data{ii,2} = dlmread(fileNames{ii});
end
number_of_experiments=length(data)/5;  %for each experiment 5 csv files are available!
for j=1:number_of_experiments
    for i=1:100
        DTs_names(j,i)=' ';
        FDs_names(j,i)=' ';
        FTs_names(j,i)=' ';
        VDs_names(j,i)=' ';
    end
end

for j=1:number_of_experiments          %making a matrix of the experiments 
    j
    w = 1+5*(j-1)
    x = 2+5*(j-1)
    y = 3+5*(j-1)
    z = 4+5*(j-1)
    
    DTs(:,j)=  data{w,2}([1:2048],2);
    DTs_ts(:,j)=  data{w,2}([1:2048],1);
    for i=1:length(data{w,1})
        if (data{w,1}(i)~='_')
            DTs_names(j,i)=data{w,1}(i);
        else
            DTs_names(j,i)='-';
        end
    end
    
    FDs(:,j)=  data{x,2}([1:2048],2);
    FDs_ds(:,j)=  data{x,2}([1:2048],1);
    for i=1:length(data{x,1})
        if (data{x,1}(i)~='_')
            FDs_names(j,i)=data{x,1}(i);
        else
            FDs_names(j,i)='-';
        end
    end
    FTs(:,j)=  data{y,2}([1:2048],2);
    FTs_ts(:,j)=  data{y,2}([1:2048],1);
    for i=1:length(data{y,1})
        
        if (data{y,1}(i)~='_')
            FTs_names(j,i)=data{y,1}(i);
        else
            FTs_names(j,i)='-';
        end
    end
    VDs(:,j)=  data{z,2}([1:2048],2);
    VDs_ds(:,j)=  data{z,2}([1:2048],1);
    for i=1:length(data{z,1})
        if (data{z,1}(i)~='_')
            VDs_names(j,i)=data{z,1}(i);
        else
            VDs_names(j,i)='-';
        end
    end
end
meane_DT = 'mean-DT'
meane_FD = 'mean-FD'
meane_FT = 'mean-FT'
meane_VD = 'mean-VD'
mean_index = number_of_experiments+1;
for i=1:length(meane_DT)
    DTs_names(mean_index,i)=meane_DT(i);
end
for i=1:length(meane_FD)
    FDs_names(mean_index,i)=meane_FD(i);
end
for i=1:length(meane_FT)
    FTs_names(mean_index,i)=meane_FT(i);
end
for i=1:length(meane_VD)
    VDs_names(mean_index,i)=meane_VD(i);
end


for i=1:2048                 %calculating the standard deviation for the corresponding points!
        STD_DTs(i) = std(DTs(i,:));
        STD_FDs(i) = std(FDs(i,:));
        STD_FTs(i) = std(FTs(i,:));
        STD_VDs(i) = std(VDs(i,:));
end

sum_Displacements= 0;
for i = 1:number_of_experiments
    
    sum_Displacements = sum_Displacements + (DTs(:,i));
    
end
mean_Displacements = sum_Displacements/number_of_experiments;
DT_std = mean(STD_DTs);
FD_std = mean(STD_FDs);
FT_std = mean(STD_FTs);
VD_std = mean(STD_VDs);
 %calculating the mean values!   

DT_mean = sum(DTs,2)/number_of_experiments;
FD_mean = sum(FDs,2)/number_of_experiments;
FT_mean = sum(FTs,2)/number_of_experiments;
VD_mean = sum(VDs,2)/number_of_experiments;

f = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(f, 'Position'); % grab original on screen size
set(f, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(f,'PaperPositionMode','auto') %set paper pos for printing
%f is for figure(1)
%plotting the force-displacement, force-time, velocity-displacement, and displacement-time diagrams for raw data.
for j=1:number_of_experiments 
    
    subplot(2,2,1)
    b(j)= plot(FDs_ds(:,j)*1000,tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000,'-','DisplayName',FDs_names(j));
    set(gca,'fontsize',20)
    xlabel('Displacement (µm)');
    ylabel('Force (mN)');
    grid on
    str_1 = sprintf('N = %d, \\sigma = %5.2f mN, N_\\chi = %d' ,number_of_experiments,FD_std*1000,FD_movingaverage);
    title(str_1);%,'interpreter','latex');
    
    
    hold on
    
    
    
    subplot(2,2,2);
    c(j)=plot(FTs_ts(:,j),tsmovavg(FTs(:,j),'t',FT_movingaverage,1)*1000,'-','DisplayName',FTs_names(j));
    set(gca,'fontsize',20)
    xlabel('Time (s)');
    ylabel('Force (mN)');
    str_2 = sprintf('N = %d, \\sigma = %5.2f mN, N_\\chi = %d',number_of_experiments,FT_std*1000,FT_movingaverage);
    title(str_2);%,'interpreter','latex');
   
    
    hold on
    grid on

    subplot(2,2,3);
    d(j)= plot(VDs_ds(:,j)*1000,tsmovavg(VDs(:,j),'t',VD_movingaverage,1),'-','DisplayName',VDs_names(j));
    set(gca,'fontsize',20)
    xlabel('Displacement (µm)');
    ylabel('Velocity (mm/s)');
    str_3 = sprintf('N = %d, \\sigma = %5.2f mm/s, N_\\chi = %d',number_of_experiments,VD_std,VD_movingaverage);
    title(str_3);%,'interpreter','latex');
    hold on
    grid on
    
    subplot(2,2,4);
    e(j)=plot(FTs_ts(:,j),tsmovavg(DTs(:,j),'t',DT_movingaverage,1)*1000,'-','DisplayName',DTs_names(j));
    set(gca,'fontsize',20)
    xlabel('Time (s)');
    ylabel('Displacement (µm)');
    str_4 = sprintf('N = %d, \\sigma = %5.2f µm, N_\\chi = %d',number_of_experiments,DT_std*1000, DT_movingaverage);
    title(str_4);%,'interpreter','latex');
    hold on
    grid on
    
end

subplot(2,2,1)
grid on
grid minor
subplot(2,2,2);
grid on
grid minor
subplot(2,2,3);
grid on
grid minor
subplot(2,2,4);
grid on
grid minor

file_name =strcat(data_path,'\','all_raws');

saveas(f, 'all_raws') % save figure
print(f,'all_raws','-dpng', '-r300')
print(f,'all_raws','-dtiff')
set(f,'Position', origSize) %set back to original dimensions
to_print_raw=strcat(file_name,'.png');

%saveFigure(figure(1),'Raw_Values');

%g is for figure(2)
g = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(g, 'Position'); % grab original on screen size
set(g, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(g,'PaperPositionMode','auto') %set paper pos for printing

%plotting the force-displacement, force-time, velocity-displacement, and displacement-time diagrams for displaced data.

for j=1:number_of_experiments 
    
    
    
    %Estimate the surface meeting force 
    val_forcy =force_threshold; %value to find
    forcies= tsmovavg(FDs(:,j),'t',FD_movingaverage,1);
    tmp1 = abs(forcies-val_forcy)
    [idxx idxx] = min(tmp1) %index of closest value
    closest_force = forcies(idxx) %closest value
    index=idxx;
    initial_distance_mm =  FDs_ds(index,j)-FDs_ds(1,j);
    initial_distance(j) = initial_distance_mm * 1000; %µm
	initial_time_s(j) = FTs_ts(index,j);
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    subplot(2,2,1)
    displacement_offset(j) = FDs_ds(1,j);
    dydx = diff([eps;tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000])./diff([eps;(FDs_ds(:,j)-displacement_offset(j))*1000-initial_distance(j)]);
    xc=1:1:2048;
    [M,I] = min(dydx(begin_diff:end_diff));
	while(abs(M)>10)
		end_diff=end_diff-100;
		[M,I] = min(dydx(begin_diff:end_diff));
	end
    I=I+begin_diff-1; %to find the real index!
    b(j)= plot((FDs_ds(:,j)-displacement_offset(j))*1000-initial_distance(j),tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000,'--','DisplayName',FDs_names(j));
    hold on;
     popo=tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000;
    %plot((FDs_ds(I,j)-displacement_offset(j))*1000-initial_distance(j),popo(I),'*','DisplayName',FDs_names(j));
    
    penetration_force(j)=max(popo(1:I));
    dropped_force (j) = min(popo(I:length(popo)));
    penetration_force_index(j)= find(popo==penetration_force(j));
    dropped_force_index(j) = find(popo==dropped_force(j));
    force_drop(j)= penetration_force(j)- dropped_force (j);
    displacement_at_penetration(j)= (FDs_ds(penetration_force_index(j),j)-displacement_offset(j))*1000-initial_distance(j) ;
	displacement_after_penetration(j)= (FDs_ds(dropped_force_index(j),j)-displacement_offset(j))*1000-initial_distance(j) ;
	drop_penetration (j)= displacement_after_penetration(j) -  displacement_at_penetration(j);
    plot((FDs_ds(penetration_force_index(j),j)-displacement_offset(j))*1000-initial_distance(j),popo(penetration_force_index(j)),'*','DisplayName',FDs_names(j));
    hold on;
     plot((FDs_ds(dropped_force_index(j),j)-displacement_offset(j))*1000-initial_distance(j),popo(dropped_force_index(j)),'*','DisplayName',FDs_names(j));
    set(gca,'fontsize',20)
    xlabel('Displacement (µm)');
    ylabel('Force (mN)');
    grid on
    % create smaller axes in top left, and plot on it
    hold on

    subplot(2,2,2);
    
    
    c(j)=plot(FTs_ts(:,j)-initial_time_s(j),(tsmovavg(FTs(:,j),'t',FT_movingaverage,1))*1000,'--','DisplayName',FTs_names(j));
    time_at_penetration(j)= FTs_ts(penetration_force_index(j),j)-initial_time_s(j);
    set(gca,'fontsize',20)
    %plot(FTs_ts([1:1900],j),(tsmovavg(FTs([1:1900],j),'t',FT_movingaverage,1))*1000,'*','DisplayName',FTs_names(j));
    %ylim([0  0.2])
    hold on

    xlabel('Time (s)');
    ylabel('Force (mN)');
    
    hold on
    grid on
    
    subplot(2,2,3);
    Vdisplacement_offset(j) = VDs_ds(1,j);
    d(j)= plot((VDs_ds(:,j)-Vdisplacement_offset(j))*1000-initial_distance(j),tsmovavg(VDs(:,j),'t',VD_movingaverage,1),'--','DisplayName',VDs_names(j));
    set(gca,'fontsize',20)
    xlabel('Displacement (µm)');
    ylabel('Velocity (mm/s)');
    hold on
    grid on
    
    subplot(2,2,4);
    Ddisplacement_offset(j) = DTs(1,j);
	         %FTs_ts(:,j)-initial_time_s(j)
    e(j)=plot(FTs_ts(:,j)-initial_time_s(j),(tsmovavg(DTs(:,j),'t',DT_movingaverage,1)-Ddisplacement_offset(j))*1000-initial_distance(j),'--','DisplayName',DTs_names(j));
    set(gca,'fontsize',20)
    xlabel('Time (s)');
    ylabel('Displacement (µm)');
    hold on
    grid on
    
end


subplot(2,2,1)

    val_forcy =force_threshold; %value to find
    forcies= tsmovavg(FD_mean,'t',FD_movingaverage,1);
    tmp1 = abs(forcies-val_forcy)
    [idxx idxx] = min(tmp1) %index of closest value
    closest_force = forcies(idxx) %closest value
    index=idxx;
    initial_distance_mm =  mean_Displacements(index)-mean_Displacements(1);
    initial_distance_mean = initial_distance_mm * 1000; %µm



b(mean_index) = plot((mean_Displacements-mean_Displacements(1))*1000-(initial_distance_mean),tsmovavg(FD_mean,'t',FD_movingaverage,1)*1000,'.-','DisplayName',FDs_names(mean_index));


%str1 = sprintf('N = %d,  \\sigma = %f, N_\\chi = %d \n $$\\bar{a}_{1}$$ = %f mN, \n Average Force Drop = %f mN \n Average Displacement at Penetration = %f µm' ,number_of_experiments, FD_std, FD_movingaverage, mean(penetration_force), mean(force_drop), mean (displacement_at_penetration));

%str40 = sprintf('$N = %d,  \\sigma = %f, N_\\chi = %d \n F_penetration_max = %f mN, F_penetration_min = %f mN \n F_drop_max = %f mN, F_drop_min = %f \n D_penetration_min = %f µm, D_penetration_max = %f µm $' ,number_of_experiments, FD_std, FD_movingaverage, ...
%        max(penetration_force), min(penetration_force), max(force_drop),  min(force_drop),  max(displacement_at_penetration),  min(displacement_at_penetration));
str50 = ['$N$ =', sprintf('%d, ', number_of_experiments),'$N_\chi$=', sprintf('%d, ', FD_movingaverage),'$\sigma$=', sprintf('%5.2f mN, ', FD_std*1000)];% , '$\sigma$=', ...
str60 = ['$\bar{F}_{penetration}$ =', sprintf('%5.2f mN, ', mean(penetration_force)), '$\bar{F}_{drop}$ =', sprintf('%5.2f mN, ', mean(force_drop))];%    , '$\bar{D}_{penetration}$ =', sprintf('%f µm,', mean (displacement_at_penetration))  ];
str70 = ['$\bar{D}_{penetration}$ =', sprintf('%5.2f ', mean(displacement_at_penetration)),'$\mu$ ' ,'$m , $'] ;

str_drop = ['$\bar{D}_{drop}$ =', sprintf('%5.2f ', mean(drop_penetration)),'$\mu$ ' ,'$m$'] ;

str_final=strvcat(str50,str60,str70,str_drop);

str80 = ['$Max_F_penetration$ =', sprintf('%5.2f mN, ', max(penetration_force)), '$Min_F_penetration$ =', sprintf('%5.2f mN, ', min(penetration_force)), '$\simga_F_penetration$ =', sprintf('%5.2f mN, ', std(penetration_force))];
str90 = ['$Max_F_drop$ =', sprintf('%5.2f mN, ', max(force_drop)), '$Min_F_drop$ =', sprintf('%5.2f mN, ', min(force_drop)), '$\simga_F_drop$ =', sprintf('%5.2f mN, ', std(force_drop))];
str100 = ['$Max_D_penetration$ =', sprintf('%5.2f microm, ', max(displacement_at_penetration)), '$Min_D_penetration$ =', sprintf('%5.2f microm, ', min(displacement_at_penetration)), '$\simga_D_penetration$ =', sprintf('%5.2f microm, ', std(displacement_at_penetration))];


ass=title(str_final,'interpreter','latex');
khar=get(get(gca,'title'),'string')
%hold on
grid on
grid minor

%ax2 = subplot(2,1,2);
%plot(x,y)
%ylim(ax2,[-10 10])

%plotting the drop of the force in a smaller box
FD_box_ax = axes('Position',[0.16 0.75 .1 .1])
box on
for j=1:number_of_experiments 
    %displacement_offset = FDs_ds(1,j);
    plot((FDs_ds(:,j)-displacement_offset(j))*1000-initial_distance(j),tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000,'--','DisplayName',FDs_names(j));
    %plot((FDs_ds(penetration_force_index(j),j)-displacement_offset(j))*1000-initial_distance(j),popo(penetration_force_index(j)),'*','DisplayName',FDs_names(j));
    %plot((FDs_ds(dropped_force_index(j),j)-displacement_offset(j))*1000-initial_distance(j),popo(dropped_force_index(j)),'*','DisplayName',FDs_names(j));
    hold on
end

%xlim ([1000 2500])
ylim (FD_box_ax,[min(penetration_force)-max(force_drop)-5 max(penetration_force)+5])
xlabel('Displacement (µm)');
ylabel('Force (mN)');
grid on 
grid minor

subplot(2,2,2);
c(mean_index)=plot(FTs_ts(:,j)-mean(initial_time_s),(tsmovavg(FT_mean,'t',FT_movingaverage,1))*1000,'.-','DisplayName',FTs_names(mean_index));

str500 = ['$N$ =', sprintf('%d, ', number_of_experiments),'$N_\chi$=', sprintf('%d, ', FT_movingaverage),'$\sigma$=', sprintf('%5.2f mN, ', FT_std*1000)];
str600 = ['$\bar{T}_{penetration}$ =', sprintf('%5.3f s, ', mean(time_at_penetration))];
str_final_time=strvcat(str500,str600);

title(str_final_time,'interpreter','latex');
str700 = ['$Max_T_penetration$ =', sprintf('%5.3f s, ', max(time_at_penetration)), '$Min_T_penetration$ =', sprintf('%5.3f s, ', min(time_at_penetration)), '$\simga_T_penetration$ =', sprintf('%5.3f s, ', std(time_at_penetration))];

grid on
grid minor

str_final_text = strvcat(str50,str60,str70,str80,str90,str100,str700);
%dlmwrite('Penetration_Characteristics.txt',str_final_text,'delimiter','');

FT_box_ax = axes('Position',[0.60 0.75 .1 .1])
box on
for j=1:number_of_experiments 
    plot((FTs_ts(:,j)-initial_time_s(j))*1000,(tsmovavg(FTs(:,j),'t',FT_movingaverage,1))*1000,'--','DisplayName',FTs_names(j));
    hold on
end

%xlim ([60 90])
ylim (FT_box_ax,[min(penetration_force)-max(force_drop)-5 max(penetration_force)+5])
xlabel('Time (ms)');
ylabel('Force (mN)');
grid on 
grid minor



subplot(2,2,3);
d(mean_index)= plot((mean_Displacements-mean(Vdisplacement_offset))*1000-mean(initial_distance),tsmovavg(VD_mean,'t',VD_movingaverage,1),'.-','DisplayName',VDs_names(mean_index));
str223 = ['$N$ =', sprintf('%d, ', number_of_experiments),'$N_\chi$=', sprintf('%d, ', VD_movingaverage),'$\sigma$=', sprintf('%5.3f mm/s, ', VD_std)];

title(str223,'interpreter','latex');
grid on
grid minor


subplot(2,2,4);
				  %FTs_ts(:,j)-mean(initial_time_s)
e(mean_index)=plot(FTs_ts(:,j)-mean(initial_time_s),(tsmovavg(DT_mean,'t',DT_movingaverage,1)-mean(Ddisplacement_offset))*1000-mean(initial_distance),'.-','DisplayName',DTs_names(mean_index));
str4 = sprintf('N = %d,  \\sigma = %5.2f µm, N_\\chi = %d \n' ,number_of_experiments, DT_std*1000, DT_movingaverage);
str224 = ['$N$ =', sprintf('%d, ', number_of_experiments),'$N_\chi$=', sprintf('%d, ', DT_movingaverage),'$\sigma$=', sprintf('%5.3f ', DT_std),'$\mu$ ' ,'$m$' ];
title(str224,'interpreter','latex');
%hold on

grid on
grid minor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


file_name =strcat(data_path,'\','all_analyzed');

    saveas(g, 'all_analyzed') % save figure
    print(g,'all_analyzed','-dpng', '-r300')
    print(g,'all_analyzed','-dtiff')
    set(g,'Position', origSize) %set back to original dimensions
to_print_analyzed=strcat(file_name,'.png');
%saveFigure(figure(2),'Displaced_Values');

%q is for figure(9)
q = figure('visible','off');
screen_size = get(0, 'ScreenSize');
origSize = get(q, 'Position'); % grab original on screen size
set(q, 'Position', [0 0 screen_size(3) screen_size(4) ] ); %set to screen size
set(q,'PaperPositionMode','auto') %set paper pos for printing


for j=1:number_of_experiments 
    xc=1:1:2048;
    displacement_offset(j) = FDs_ds(1,j);
    dydx = diff([eps;tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000])./diff([eps;(FDs_ds(:,j)-displacement_offset(j))*1000-initial_distance(j)]);
	set(gca,'fontsize',20)
    subplot(2,1,2);
    plot(xc,tsmovavg(FDs(:,j),'t',FD_movingaverage,1)*1000,'Linewidth', 2)
	hold on
	plot(penetration_force_index(j),popo(penetration_force_index(j)),'o','DisplayName',FDs_names(j),'Markersize',10,'Linewidth', 3);
    hold on;
    plot(dropped_force_index(j),popo(dropped_force_index(j)),'o','DisplayName',FDs_names(j),'Markersize',10,'Linewidth', 3);
    hold on
    grid on
    grid minor
    xlabel('index');
    ylabel('Force (mN)');
	set(gca,'fontsize',30)
    subplot(2,1,1);
    plot(xc,tsmovavg(dydx,'t',FD_movingaverage,1),'Linewidth', 2)
    xlabel('index');
    ylabel('derivative value');
	set(gca,'fontsize',30)
    grid on 
    grid minor
    hold on
end
grid on
grid minor
 subplot(2,1,1);
grid on
grid minor
 subplot(2,1,2);
grid on
grid minor
 subplot(2,1,1);
grid on
grid minor
file_name =strcat(data_path,'\','derivative_index');
saveas(q, 'derivative_index') % save figure
print(q,'derivative_index','-dpng', '-r300')
print(q,'derivative_index','-dtiff')
set(q,'Position', origSize) %set back to original dimensions
to_print_derivative=strcat(file_name,'.png');

%openfig('derivative_index.fig','visible');


read_filename='readme.txt';
the_path =strcat(read_filename);
dd = dir(read_filename);                           
fileName = {dd.name}; 
file    = memmapfile( char(fileName), 'writable', true );
comma   = uint8(',');
point   = uint8('.');
file.Data( transpose( file.Data==comma) ) = point;
[parameters, values] = textread(read_filename,'%s %s', 11);

IndexC = strfind(parameters, 'Diameter');
Index = find(not(cellfun('isempty', IndexC)));
diameter = str2double(values(Index));

IndexC = strfind(parameters, 'Number');
Index = find(not(cellfun('isempty', IndexC)));
number_of_microneedles = str2double(values(Index));

IndexC = strfind(parameters, 'mm/s');
Index = find(not(cellfun('isempty', IndexC)));
Velocity = str2double(values(Index));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5*
number_of_experiments

mean(penetration_force)
std(penetration_force)
max(penetration_force) 
min(penetration_force)


mean(force_drop)
std(force_drop)
max(force_drop)
min(force_drop) 


mean(displacement_at_penetration)
std(displacement_at_penetration)
max(displacement_at_penetration)
min(displacement_at_penetration)


mean(time_at_penetration)
std(time_at_penetration)
max(time_at_penetration)
min(time_at_penetration)

cd('C:/')
existence_fin=exist('final_analysis','dir')
if(existence_fin==0)
	mkdir('final_analysis')
	end
%%%%%%%%%%%%%%%%%%%%%The average and standard deviation values along with the microneedle properties are saved in a csv file named final_analysis with the velocity and number of microneedles which is further used for bar graph production%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
final_analysis_filename= strcat('final_analysis_',num2str(number_of_microneedles),'mn_','V_',num2str(Velocity),'mms-1','.csv');
final_analysis_path='C:\final_analysis\';
final_analysis_exist=strcat(final_analysis_path,final_analysis_filename);
a=exist(final_analysis_exist, 'file')
cd(final_analysis_path)
if(a==0)
	fid = fopen(final_analysis_filename, 'w');
	fprintf(fid, 'Tip Diameter; Penetration Force Average (mN); Penetration Force Standard Deviation (mN) ; Penetration Force Maximum (mN); Penetration Force Minimum (mN);Force Drop Average (mN); Force Drop Standard Deviation (mN) ; Force Drop Maximum (mN) ; Force Drop Minimum (mN); Displacement at Penetration Average (µm); Displacement at Penetration Standard Deviation (µm); Displacement at Penetration Maximum (µm); Displacement at Penetration Minimum (µm) ;Time at Penetration Average (s); Time at Penetration Standard Deviation (s) ; Time at Penetration Maximum (s) ; Time at Penetration Minimum (s)); Number of Experiments');
	fclose(fid)
end

fid = fopen(final_analysis_filename, 'a');
fprintf(fid,strcat('\n',strrep(num2str(diameter),'.', ','),';',strrep(num2str(mean(penetration_force)),'.', ','),';', strrep(num2str(std(penetration_force)),'.', ',') ,';', strrep(num2str(max(penetration_force)),'.', ',') ,';', strrep(num2str(min(penetration_force)) ,'.', ','),';', strrep(num2str(mean(force_drop)) ,'.', ','),';', strrep(num2str(std(force_drop)) ,'.', ','),';', strrep(num2str(max(force_drop)),'.', ',') ,';', strrep(num2str(min(force_drop)),'.', ','),';', strrep(num2str(mean(displacement_at_penetration)),'.', ','),';', strrep(num2str(std(displacement_at_penetration)),'.', ','),';', strrep(num2str(max(displacement_at_penetration)),'.', ','),';', strrep(num2str(min(displacement_at_penetration)),'.', ','),';', strrep(num2str(mean(time_at_penetration)) ,'.', ','),';', strrep(num2str(std(time_at_penetration)) ,'.', ','),';', strrep(num2str(max(time_at_penetration)),'.', ',') ,';',strrep(num2str(min(time_at_penetration)),'.', ','),';',strrep(num2str(number_of_experiments),'.', ',')));
fclose(fid);

