--[[
 SIMION Lua workbench program for octupole simulation.
 This oscillates octupole rod potentials
 (and also updates PE display periodically).

 D.Manura-2007-09
 (c) 2007 Scientific Instrument Services, Inc. (Licensed under SIMION 8.0)
--]]

simion.workbench_program()

-- Variables adjustable during flight:

adjustable pe_update_each_usec      = 0.05   -- potential energy display
                                             -- update period (microsec)
                                             -- (for display purposes only)

-- Variables adjustable only at beginning of flight:

adjustable effective_radius_mm      = 8.33    -- half the minimum distance between
                                             -- opposite rods (mm)
adjustable phase_angle_deg          = 0.0    -- entry phase angle of ion (deg)
adjustable _frequency_hz            = 1.2E6  -- RF frequency (Hz)
adjustable _a4
adjustable _q4


function compute_voltages()

  local q = ion_charge * 1.602176462*10^-19
  local m = ion_mass * 1.66053885*10^-27

  local omega = _frequency_hz * (1E-6 * 2 * math.pi)
  local theta = phase_angle_deg * (math.pi / 180)


  local r0 = effective_radius_mm / 1000 -- (m/mm)
  adjustable _dcvolts = _a4 * m * omega^2 * r0^2 / (8 * q)      -- DC voltage for octupole;
                                                           --    typically zero for RF-only octupoles
  adjustable _rfvolts = _q4 * m * omega^2 * r0^2 / (4 * q)           -- RF voltage for octupole
  return _dcvolts, _rfvolts
end

-- internal variables
local last_pe_update = 0.0 -- last potential energy surface update time (usec)


function segment.flym()
  -- Simulate all parameterizations.
  -- Results are summarized in results.csv.
  -- enable 'output' directory exists
  os.execute("mkdir output")

  -- cleanup any old files.
  os.remove("output/q_a_optimization.out")

  -- Open output files.
  local results_file = assert(io.open("output/q_a_optimization.out", "w")) -- write mode
  
  run_number = 0
  -- loop over all (a, q) values
  for ai= 0.220, 0.250, 0.010
  do
    for qi = 0.700, 0.720, 0.010 
    do
      a4 = ai
      q4 = qi

      -- run for each (ai, qi) value
      
      run()
      run_number = run_number + 1

      -- retrieve transmission efficiency for all mass_numbers 
      num_ions, num_hits, transmission = segment.terminate_run()
      results_file:write("run_number, a, q \n")
      results_file:write(run_number.. ", " .. a4 .. ", " .. q4 .." \n\n")
      results_file:write("mass_number, num_ions, fraction_ions, num_hits, efficiency (%)\n")

      results_file:write("total ions flown: ".. total_ions .."\n")
    
      for i=1,7 do
        print(transmission[i])
        results_file:write(mass_number[i].. ", " .. num_ions[i].. ", " .. num_ions[i]/total_ions .. ", " .. num_hits[i] .. ", " .. transmission[i] .."\n")
        results_file:flush()
      end

      results_file:write("---------------------------------------- \n\n")
    end
  end
end


-- SIMION segment called by SIMION to set adjustable electrode voltages
-- in the current potential array instance.
-- NOTE: this is called frequently, multiple times per time-step (by
-- Runge-Kutta), so performance concerns here can be important.
function segment.fast_adjust()

_dcvolts, _rfvolts = compute_voltages()

  local tempvolts =
    sin(ion_time_of_flight * omega + theta) * _rfvolts + _dcvolts

  -- Apply adjustable voltages to rod electrodes.
  adj_elect04 =   tempvolts
  adj_elect06 =   tempvolts
  adj_elect05 = - tempvolts
  adj_elect07 = - tempvolts
end
-- See also the README.html for how memory usage might be further
-- reduced by 1/3 or 2/3 by sharing electrode solution arrays.


-- SIMION segment called by SIMION after every time-step.
local is_initialized
function segment.other_actions()
  if not is_initialized then
    -- Convert to SI units.
    local q = ion_charge * 1.602176462*10^-19 -- (C/e)
    local m = ion_mass * 1.66053886*10^-27    -- (kg/u)
       
    -- Print stability constants.
    print(string.format("m/z=%g,a4=%g,q4=%g", ion_mass/ion_charge, a4, q4))
    is_initialized = true -- only execute once
  end

  -- Update potential energy surface display periodically.
  -- The performance overhead of this in non-PE views is only a few percent.
  -- NOTE: the value inside abs(...) can be negative when a new ion is flown.
  if abs(ion_time_of_flight - last_pe_update) >= pe_update_each_usec then
    last_pe_update = ion_time_of_flight
    sim_update_pe_surface = 1    -- Request a PE surface display update.
  end
end

-- SIMION segment called by SIMION to override time-step size on each time-step.
function segment.tstep_adjust()
   -- Keep time step size <= X usec.
   ion_time_step = min(ion_time_step, 0.1)  -- X usec
end



--- Added a posteriori from count.lua

-- Counters for flown ions with mass_numbers.

local total_ions = 0
local num_ions = {0, 0, 0, 0, 0, 0, 0, 0}
local num_hits = {0, 0, 0, 0, 0, 0, 0, 0}
local ystart = {}
local mass_number = {133, 134, 135, 136, 137, 138, 139, 140, 141}

-- called on start of each run...
function segment.initialize_run()
  -- Reset the counter before each rerun (only needed if Rerun is enabled).
  --num_hits = 0
  segment.fast_adjust()
 end

-- called on each particle initialization inside a PA instance...
function segment.initialize()
  -- Infer total number of particles flown.  [*1]
  --- num_ions = ion_number  this way particles not initialize are accounted for as well, thus underestimating the transmission efficiency
  
    -- Optionally store data on this particle's starting conditions.  [*2]
    ystart[ion_number] = ion_py_mm
end


  -- called on each particle termination inside a PA instance...
function segment.terminate()
    -- Print data on each splat.
  --print('splat:', 'y_begin=', ystart[ion_number], 'y_end=', ion_py_mm)
    for i=1,8,1 do
    
    if ion_mass >= mass_number[i] and ion_mass <= mass_number[i+1] then
      num_ions[i] = num_ions[i] + 1   -- this way only the particles that splat are accounted for
      --print(i, ",", mass_number[i], "ion_pz_mm = ", ion_pz_mm, "num_hits[i] = ")
    
      if ion_pz_mm > 220 and ion_px_mm < 38 and ion_py_mm < 38 then
        num_hits[i] = num_hits[i] + 1
      end
    end
  end  
  total_ions = total_ions + 1
end


-- called on end of each run...
function segment.terminate_run()
  -- Print summary data at end of run.
  local transmission = {}
  for i=1,8,1 do
    if num_ions[i] > 0 then
      transmission[i] = 100 * num_hits[i] / num_ions[i]
  else
    transmission[i] = 0
  end
    print('num_ions=',  num_ions[i])
    print('num_hits=',       num_hits[i])
    print('efficiency (%)=', transmission[i])
  end
  return num_ions, num_hits, transmission
end


