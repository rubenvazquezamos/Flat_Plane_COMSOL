function [pressure] = QR_5(Freq,Probe,File)

%-------------------------------------------------------------------------%
%% CALL COMSOL
disp(' -- Call COMSOL')
%
% QR_5_toMATLAB.m
%
% Model exported on Aug 20 2025, 15:22 by COMSOL 6.0.0.318.

import com.comsol.model.*
import com.comsol.model.util.*
model = ModelUtil.create('Model');

model.component.create('comp1', true); %Create component 1
%-------------------------------------------------------------------------%
%% PARAMETERS
disp(' -- Sending parameters')
%-------------------------------------------------------------------------%

model.param.set('L', '0.52[m]', 'Width of the diffuser');
model.param.set('Lw', '0.09[m]', 'Width of one well');
model.param.set('Li', '0.01[m]', 'Width in between wells');
model.param.set('H', '0.3[m]', 'Height of the diffuser');
model.param.set('d1', '0[m]', 'Depth of well 1');
model.param.set('d2', '0.068[m]', 'Depth of well 2');
model.param.set('d3', '0.272[m]', 'Depth of well 3');
model.param.set('d4', '0.272[m]', 'Depth of well 4');
model.param.set('d5', '0.068[m]', 'Depth of well 5');
model.param.set('d6', '0.000001[m]', 'Depth of well 6');
model.param.set('r_air', '(L/2)*20', 'Radius of the air domain (for single diffuser)');
model.param.set('r0', '100[m]', 'Evaluation distance');
model.param.set('Hair', '1[m]', 'Height of the air domain');
model.param.set('Hpml', '0.2[m]', 'Thickness of the PML');
model.param.group.create('par2');
model.param('par2').set('c0', '343[m/s]', 'Speed of sound');
model.param('par2').set('rho0', '1.225[kg/m^3]', 'Density');
model.param('par2').set('f_min', num2str(Freq.f_min));
model.param('par2').set('f_max', num2str(Freq.f_max));
model.param('par2').set('df', num2str(Freq.df));
model.param('par2').set('phi', '0 [rad]', 'phase angle');
model.param.label('Geometry');
model.param('par2').label('Physics');

%-------------------------------------------------------------------------%
%% COMSOL PHYSICAL CONSTANTS
disp(' -- Set global constants')
%-------------------------------------------------------------------------%

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').propertyGroup('def').func.create('eta', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('Cp', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('rho', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('k', 'Piecewise');
model.component('comp1').material('mat1').propertyGroup('def').func.create('cs', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an1', 'Analytic');
model.component('comp1').material('mat1').propertyGroup('def').func.create('an2', 'Analytic');
model.component('comp1').material('mat1').propertyGroup.create('RefractiveIndex', 'Refractive index');
model.component('comp1').material('mat1').propertyGroup.create('NonlinearModel', 'Nonlinear model');
model.component('comp1').material('mat1').propertyGroup.create('idealGas', 'Ideal gas');
model.component('comp1').material('mat1').propertyGroup('idealGas').func.create('Cp', 'Piecewise');
%-------------------------------------------------------------------------%
%% VARIABLES
disp(' -- Sending variables')
%-------------------------------------------------------------------------%
model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('p_inc', '1[Pa]*exp(i*2*pi*freq/c0*(sin(theta0)*x-cos(theta0)*y))', 'Incident plane wave');
model.component('comp1').variable('var1').set('p_inf', '1[Pa]*exp(i*2*pi*freq/c0*(sin(theta0)*x+cos(theta0)*y))', 'Reflected plane wave from infinite baffle');
model.component('comp1').variable('var1').set('p_scat', 'p_inf+acpr.p_s', 'Scattered pressure');
model.component('comp1').variable('var1').set('p_scat_ext', 'pext(x,y)', 'Scattered pressure in the exterior field');

%-------------------------------------------------------------------------%
%% GEOMETRY
disp(' -- Build geometry')
%-------------------------------------------------------------------------%

model.component('comp1').geom.create('geom1', 2);
model.component('comp1').mesh.create('mesh1');
model.component('comp1').geom('geom1').create('c1', 'Circle');
model.component('comp1').geom('geom1').feature('c1').set('r', 'r_air');
model.component('comp1').geom('geom1').feature('c1').set('angle', 180);
model.component('comp1').geom('geom1').create('r1', 'Rectangle');
model.component('comp1').geom('geom1').feature('r1').active(false);
model.component('comp1').geom('geom1').feature('r1').set('pos', {'-L/2+Li' '-d1'});
model.component('comp1').geom('geom1').feature('r1').set('size', {'Lw' 'd1'});
model.component('comp1').geom('geom1').create('r2', 'Rectangle');
model.component('comp1').geom('geom1').feature('r2').set('pos', {'-L/2+2*Li+Lw' '-d2'});
model.component('comp1').geom('geom1').feature('r2').set('size', {'Lw' 'd2'});
model.component('comp1').geom('geom1').create('r3', 'Rectangle');
model.component('comp1').geom('geom1').feature('r3').set('pos', {'-L/2+3*Li+2*Lw' '-d3'});
model.component('comp1').geom('geom1').feature('r3').set('size', {'Lw' 'd3'});
model.component('comp1').geom('geom1').create('r4', 'Rectangle');
model.component('comp1').geom('geom1').feature('r4').set('pos', {'-L/2+4*Li+3*Lw' '-d4'});
model.component('comp1').geom('geom1').feature('r4').set('size', {'Lw' 'd4'});
model.component('comp1').geom('geom1').create('r5', 'Rectangle');
model.component('comp1').geom('geom1').feature('r5').set('pos', {'-L/2+5*Li+4*Lw' '-d5'});
model.component('comp1').geom('geom1').feature('r5').set('size', {'Lw' 'd5'});
model.component('comp1').geom('geom1').create('uni1', 'Union');
model.component('comp1').geom('geom1').feature('uni1').set('intbnd', false);
model.component('comp1').geom('geom1').feature('uni1').selection('input').set({'c1' 'r1' 'r2' 'r3' 'r4' 'r5'});
model.component('comp1').geom('geom1').create('pt1', 'Point');
model.component('comp1').geom('geom1').feature('pt1').set('p', {'-L/2+Li' '0'});
model.component('comp1').geom('geom1').nodeGroup.create('grp1');
model.component('comp1').geom('geom1').nodeGroup('grp1').label('Wells');
model.component('comp1').geom('geom1').nodeGroup('grp1').placeAfter('c1');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r1');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r2');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r3');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r4');
model.component('comp1').geom('geom1').nodeGroup('grp1').add('r5');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

flat = questdlg("do you want to simulate a flat plane?");
switch flat
    case 'Yes'
    model.component('comp1').geom('geom1').feature('r2').active(false);
    model.component('comp1').geom('geom1').feature('r3').active(false);
    model.component('comp1').geom('geom1').feature('r4').active(false);
    model.component('comp1').geom('geom1').feature('r5').active(false);
    % model.component('comp1').geom('geom1').feature('r6').active(false);

    
    case 'No'
end

%%% EXTRA CODE: visualise the geometry before going further with
% approval input with the [Enter] key.
fig = figure();
clf(fig)
set(fig,'Renderer','opengl');
mphgeom(model,'geom1','facealpha',0.5);
box on;
kk = input('Is the geometry valid?');clear kk;
%-------------------------------------------------------------------------%

%Define explicit selection for far field coordinates
model.component('comp1').selection.create('sel1', 'Explicit');
model.component('comp1').selection('sel1').geom('geom1', 1);
model.component('comp1').selection('sel1').set([19 20]);
model.component('comp1').selection('sel1').label('Far Field');

%-------------------------------------------------------------------------%
%% PHYSICS
disp(' -- Implementing physics')
%-------------------------------------------------------------------------%

model.component('comp1').physics.create('acpr', 'PressureAcoustics', 'geom1');

model.component('comp1').physics('acpr').create('bpf1', 'BackgroundPressureField', 2);
model.component('comp1').physics('acpr').feature('bpf1').selection.set([1]);
model.component('comp1').physics('acpr').create('efc1', 'ExteriorFieldCalculation', 1);
model.component('comp1').physics('acpr').feature('efc1').selection.set([19 20]);
model.component('comp1').physics('acpr').feature('bpf1').set('p', 'p_inc');
model.component('comp1').physics('acpr').feature('bpf1').set('dir', [0; -1; 0]);
model.component('comp1').physics('acpr').feature('bpf1').set('phi', 'phi');
model.component('comp1').physics('acpr').feature('bpf1').set('pamp', 1);
model.component('comp1').physics('acpr').feature('bpf1').set('c_mat', 'from_mat');

%% PML
disp(' -- Set PML')
%-------------------------------------------------------------------------%
model.component('comp1').physics('acpr').create('pmb1', 'PerfectlyMatchedBoundary', 1);
model.component('comp1').physics('acpr').feature('pmb1').selection.set([19 20]);
model.component('comp1').physics('acpr').create('pmb2', 'PerfectlyMatchedBoundary', 1);
model.component('comp1').physics('acpr').feature('pmb2').selection.set([1 18]);
model.component('comp1').physics('acpr').feature('pmb2').set('directionType', 'normal');

%-------------------------------------------------------------------------%
%% MATERIAL
disp(' -- Define materials')
%-------------------------------------------------------------------------%

model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').set('family', 'air');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('pieces', {'200.0' '1600.0' '-8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('eta').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('expr', 'pA*0.02897/R_const[K*mol/J]/T');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('fununit', 'kg/m^3');
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('argunit', {'Pa' 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('rho').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '293.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('pieces', {'200.0' '1600.0' '-0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4'});
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('def').func('k').set('fununit', 'W/(m*K)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('expr', 'sqrt(1.4*R_const[K*mol/J]/0.02897*T)');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('fununit', 'm/s');
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('cs').set('plotargs', {'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('funcname', 'alpha_p');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('expr', '-1/rho(pA,T)*d(rho(pA,T),T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('args', {'pA' 'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('fununit', '1/K');
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('argunit', {'Pa' 'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an1').set('plotargs', {'pA' '101325' '101325'; 'T' '273.15' '373.15'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('funcname', 'muB');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('expr', '0.6*eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('args', {'T'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('fununit', 'Pa*s');
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('argunit', {'K'});
model.component('comp1').material('mat1').propertyGroup('def').func('an2').set('plotargs', {'T' '200' '1600'});
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', '');
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', '');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)' '0' '0' '0' 'alpha_p(pA,T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('molarmass', '0.02897[kg/mol]');
model.component('comp1').material('mat1').propertyGroup('def').set('bulkviscosity', 'muB(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('dynamicviscosity', 'eta(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('def').set('electricconductivity', {'0[S/m]' '0' '0' '0' '0[S/m]' '0' '0' '0' '0[S/m]'});
model.component('comp1').material('mat1').propertyGroup('def').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('def').set('density', 'rho(pA,T)');
model.component('comp1').material('mat1').propertyGroup('def').set('thermalconductivity', {'k(T)' '0' '0' '0' 'k(T)' '0' '0' '0' 'k(T)'});
model.component('comp1').material('mat1').propertyGroup('def').set('soundspeed', 'cs(T)');
model.component('comp1').material('mat1').propertyGroup('def').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('def').addInput('pressure');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', '');
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('n', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('RefractiveIndex').set('ki', {'0' '0' '0' '0' '0' '0' '0' '0' '0'});
model.component('comp1').material('mat1').propertyGroup('NonlinearModel').set('BA', '(def.gamma+1)/2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').label('Piecewise 2');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('arg', 'T');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('pieces', {'200.0' '1600.0' '1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4'});
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('argunit', 'K');
model.component('comp1').material('mat1').propertyGroup('idealGas').func('Cp').set('fununit', 'J/(kg*K)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('Rs', 'R_const/Mn');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('heatcapacity', 'Cp(T)');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('ratioofspecificheat', '1.4');
model.component('comp1').material('mat1').propertyGroup('idealGas').set('molarmass', '0.02897');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('temperature');
model.component('comp1').material('mat1').propertyGroup('idealGas').addInput('pressure');
model.component('comp1').material('mat1').materialType('nonSolid');

%-------------------------------------------------------------------------%
%% STUDY
disp(' -- Creating study')
%-------------------------------------------------------------------------%

model.sol.create('sol1');
model.study.create('std1');
model.study('std1').create('freq', 'Frequency');
model.study('std1').feature('freq').set('useadvanceddisable', true);
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('p1', 'Parametric');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');
model.study('std1').feature('freq').set('plist', 'range(f_min,df,f_max)');
model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').label('Compile Equations: Frequency Domain');
model.sol('sol1').feature('v1').label('Dependent Variables 1.1');
model.sol('sol1').feature('v1').set('clistctrl', {'p1'});
model.sol('sol1').feature('v1').set('cname', {'freq'});
model.sol('sol1').feature('v1').set('clist', {'range(f_min,df,f_max)'});
model.sol('sol1').feature('s1').label('Stationary Solver 1.1');
model.sol('sol1').feature('s1').feature('dDef').label('Direct 1');
model.sol('sol1').feature('s1').feature('aDef').label('Advanced 1');
model.sol('sol1').feature('s1').feature('aDef').set('complexfun', true);
model.sol('sol1').feature('s1').feature('p1').label('Parametric 1.1');
model.sol('sol1').feature('s1').feature('p1').set('pname', {'freq'});
model.sol('sol1').feature('s1').feature('p1').set('plistarr', {'range(f_min,df,f_max)'});
model.sol('sol1').feature('s1').feature('p1').set('punit', {'Hz'});
model.sol('sol1').feature('s1').feature('p1').set('pcontinuationmode', 'no');
model.sol('sol1').feature('s1').feature('p1').set('preusesol', 'auto');
model.sol('sol1').feature('s1').feature('p1').set('uselsqdata', false);
model.sol('sol1').feature('s1').feature('fc1').label('Fully Coupled 1.1');


%% MESH
disp(' -- Build mesh')
%-------------------------------------------------------------------------%

model.component('comp1').mesh('mesh1').create('size1', 'Size');
model.component('comp1').mesh('mesh1').create('se1', 'SizeExpression');
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').create('bl1', 'BndLayer');
model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('size1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('se1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('se1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('bl1').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('bl1').selection.set([1]);
model.component('comp1').mesh('mesh1').feature('bl1').create('blp1', 'BndLayerProp');
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp1').selection.set([19 20]);

model.component('comp1').mesh('mesh1').feature('size1').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size1').set('hmax', 'c0/f_max/5');
model.component('comp1').mesh('mesh1').feature('size1').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size1').set('hmin', 'c0/f_max/10');
model.component('comp1').mesh('mesh1').feature('size1').set('hminactive', false);

model.component('comp1').mesh('mesh1').feature('se1').set('evaltype', 'initialexpression');
model.component('comp1').mesh('mesh1').feature('se1').set('sizeexpr', 'subst(real(acpr.c_c),acpr.freq,freqmax)/freqmax/5');

model.component('comp1').mesh('mesh1').feature('bl1').set('sharpcorners', 'trim');
model.component('comp1').mesh('mesh1').feature('bl1').set('smoothtransition', false);
model.component('comp1').mesh('mesh1').feature('bl1').feature('blp1').set('blnlayers', 1);

model.component('comp1').mesh('mesh1').feature('se1').set('studystep', 'std1/freq');
%-------------------------------------------------------------------------%

%% SOLVE
disp(' -- Solving')

%%% EXTRA CODE: opens a window to show the solver's progress
ModelUtil.showProgress(true);
% runs the solver
model.sol('sol1').runAll;

%% PLOT
disp(' -- Generate plots')
%-------------------------------------------------------------------------%
model.result.create('pg3', 'PolarGroup');
model.result('pg3').create('rp1', 'RadiationPattern');
model.result('pg3').feature('rp1').set('data', 'dset1');
model.result('pg3').feature('rp1').set('expr', 'p_scat_ext');
model.result('pg3').set('zeroangle', 'up');
model.result('pg3').feature('rp1').set('phidisc', 180);
model.result('pg3').feature('rp1').set('anglerestr', 'manual');
model.result('pg3').feature('rp1').set('phimin', -90);
model.result('pg3').feature('rp1').set('phirange', 180);
model.result('pg3').feature('rp1').set('circle', 'manual');
model.result('pg3').feature('rp1').set('radius', 'r0');
model.result('pg3').feature('rp1').set('refdir2', [0 1]);
%-------------------------------------------------------------------------%
%% GET DATA/PROBING
disp(' -- Probing')
%-------------------------------------------------------------------------%
pressure = mphinterp(model,'acpr.p_s','coord',Probe.Coordinates);

%% SAVE
disp(' -- Save model')
%-------------------------------------------------------------------------%
% 
model.sol('sol1').clearSolutionData; %Clear solution data
model.mesh.clearMeshes; %Clear meshes
mphsave(model,[File.Path,filesep,File.Tag,File.Extension]);

disp(' -- DONE')

end