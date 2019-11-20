simion.workbench_program()

--listing all local and adjustable variables which here
 
adjustable _percent_tune          =    97.0  		-- percent of optimum tune.
                                             		-- (typically just under 100)
adjustable mass_in_u   =   134.0  			-- mass (u)
local mass = mass_in_u * 1.66054E-27                   	-- mass(kg)
	
adjustable pe_update_each_usec      = 5 		-- potential energy display update (usec)                                    

adjustable radius_in_mm   = 8.331   			--effective radius in mm
                                             
adjustable phase_angle_deg          = 0.0    		-- quad entry phase angle of ion (deg)
adjustable frequency_hz             = 1.2E6  		-- RF frequency of quad (Hz)
adjustable charge= 1					--charge state of the ion
adjustable q = charge * 1.602176462E-19        		--ion charge in Coulomb, 1e for testing purposes


local omega 	= frequency_hz * (2 * math.pi)   -- frequency_hz (reexpressed in units of radians/usec)
local theta 	= phase_angle_deg * (math.pi / 180)     -- phase_angle_deg (reexpressed in units of radians)
local last_pe_update = 0.0 				-- last potential energy surface update time (usec)
local _a4 						-- a value from stability diag
local _q4						-- q value from stability diag
local r0=radius_in_mm/1000				--effective radius of QMF rods in meter



local zstart = {}					--zstart
local mass_number = {133, 134, 135, 136, 137, 138, 139, 140, 141} -- list of mass numbers to investigate 
local run_number = 1 				--counter for changing the mass after each a4 and q4 loops are finished
local num_particles					--number of ions
local transmission = 0				--initialize transmission

--FLY'M SEGMENT-----

function segment.flym()

  -- Results are summarized in results.csv.
  -- enable 'output' directory exists
  os.execute("mkdir output")

  -- delete any old files.
  os.remove("output/q_a_optimization.csv")

  -- Open output files.
  local results_file = assert(io.open("output/q_a_optimization.out", "w")) -- write mode
  	
  -- loop over all masses and (a, q) values (one might want to add an automatic adjustment in case the list length changes)
for i=1,8,1 do
  for ai= 0.220, 0.250, 0.01
  do
    for qi = 0.700, 0.720, 0.01 
    do
      _a4 = ai
      _q4 = qi

      -- run for each (ai, qi) value, mass is changed in initialize segment
      
      run()
      
		
	num_ions, num_hits, transmission = segment.terminate_run()
	--write stuff to file
      results_file:write("run_number, a, q \n")
      results_file:write(run_number.. ", " .. _a4 .. ", " .. _q4 .." \n\n")
      

     results_file:write("mass_number, num_ions, num_hits, transmission (%)\n")
    
     
       results_file:write(mass_number[run_number+1].. ", " .. num_particles.. ", "  .. num_hits .. ", " .. transmission .."\n")
       results_file:flush()

      results_file:write("---------------------------------------- \n\n")
		end
	end
    run_number = run_number + 1 --counter+1 to change mass (down in initialize segment)
	end
end



-- FAST ADJUST SEGMENT

function segment.fast_adjust()
    --calculate and apply voltages, can probably replaced with its own function again, for testing purposes I put the calculation here

    local rfvolts = _q4 * mass * omega^2 * r0^2 / (4 * q)
    local dcvolts = _a4 * mass * omega^2 * r0^2 / (8 * q) 
    local tempvolts = sin(ion_time_of_flight *1E-6 *  omega+ theta) * rfvolts + dcvolts

    adj_elect07 =tempvolts 
    adj_elect05 =tempvolts
    adj_elect04 =-tempvolts
    adj_elect06 =-tempvolts
    
	end


--OTHER ACTIONS SEGMENT

-- SIMION segment called by SIMION after every time-step.
function segment.other_actions()
 --the PE update etc might not work, did not check
  	if abs(ion_time_of_flight - last_pe_update) >= pe_update_each_usec then
    		last_pe_update = ion_time_of_flight
    		sim_update_pe_surface = 1  -- Request a PE surface display update.
  	end
	
end


-- INITIALIZE AND INITIALIZE RUN SEGMENTS

function segment.initialize()
--called each time a particle is initialized inside a PA instance 

print('current mass=',mass_number[run_number]) --for testing, wanted to see mass number for each new a,q set. Can be disabled

num_particles = ion_number --create non system internal variable equal to number of ions in simulation
zstart[ion_number] = ion_pz_mm --Z starting point, adjust according to your geometry/ direction of ion movement

ion_mass= mass_number[run_number+1] --change ion_mass with counter from above

 
end


function segment.initialize_run()
--called exactly once prior to the start of each run 

  num_hits = 0
 	print('running with q4=',_q4) --also for testing, to see things are working
	print('running with a4=',_a4) --also for testing, to see things are working
 end


-- TERMINATE AND TERMINATE RUN SEGMENTS

function segment.terminate() 
    
 -- print('splat:', 'z_begin=', zstart[ion_number], 'z_end=', ion_pz_mm) -- Print data on each splat (for tests), can be ignored.

  -- Count particles that splat within some region of volume
  if ion_px_mm <= 38 and ion_py_mm < 38 and ion_pz_mm > 220 then --count transmitted ions. Be careful that your ions are created within the PA (workbench might be bigger). In my case, the PA(=QMF) ends at z=0mm, so particles approach max 1E-7mm 
    
      num_hits = num_hits + 1
    
  end
end


function segment.terminate_run()
 -- Print summary data at end of run and calculate transmission.
  local transmission= 100 * num_hits / num_particles
  print('num_particles=',  num_particles)
  print('num_hits=',       num_hits)
  print('efficiency (%)=', transmission)
  return num_ions, num_hits, transmission
end
