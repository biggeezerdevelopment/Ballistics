// a Ballistic library for Go

package ballistic

import (
	"encoding/csv"
	"fmt"
	"math"
	"os"
)

// MAXRANGE is the MAXIMUM range that can be calculated
const MAXRANGE = 50001

// GRAVITY - Get standard gravity and convert from m/s^2 to ft/s^2
const GRAVITY = -float64(9.80665) * 3.28084

// EARTH_ROTATION - Speed of earth's rotation (rad/s)
const EARTH_ROTATION = 0.00007292

// DragFunction - Different drag functions that can be used
type DragFunction int

// Different available drag functions
const (
	G1 DragFunction = iota + 1
	G2
	G3
	G4
	G5
	G6
	G7
	G8
)

// Solution is what we return from a calculation
type Solution struct {
	Range      int     `json:"range"`       // Range in yards
	RawRange   float64 `json:"raw_range"`   // Range in yards for calculations
	Drop       float64 `json:"drop"`        // Bullet drop in inches
	DropMOA    float64 `json:"drop_moa"`    // Bullet drop in MOA
	Time       float64 `json:"time"`        // Time of flight in seconds
	Windage    float64 `json:"windage"`     // Wind compensation in inches
	WindageMOA float64 `json:"windage_moa"` // Wind compensation in inches
	Energy     float64 `json:"energy"`      // Kinetic Energy in ft*lbs
	Velocity   float64 `json:"velocity"`    // Velocity in fps
	VelocityX  float64 `json:"velocity_x"`  // Velocity X-component
	VelocityY  float64 `json:"velocity_y"`  // Velocity Y-component
}

// Very verbose String function, mostly for debugging
func (s *Solution) String() string {
	return fmt.Sprintf("Range(yrd):%d RawRange(yrd):%0.2f Drop(in):%0.2f Drop(MOA):%0.2f TimeOfFlight(s):%0.2f Windage(in):%0.2f Windage(MOA):%0.2f Velocity(ft/s):%0.2f VelocityX(ft/s):%0.2f VelocityY(ft/s):%0.2f", s.Range, s.RawRange, s.Drop, s.DropMOA, s.Time, s.Windage, s.WindageMOA, s.Velocity, s.VelocityX, s.VelocityY)
}

// Bullet contains all the needed characteristics of a bullet
type Bullet struct {
	Caliber float64      `json:"caliber"` // Diameter of bullet in inches
	Weight  float64      `json:"weight"`  // Weight of bullet in grains
	BC      *Coefficient `json:"bc"`      // Ballistic Coefficient of bullet
	Length  float64      `json:"length"`  // bullet length in inches
}

// Coefficient defines the value and type of drag function
type Coefficient struct {
	Value    float64      `json:"value"`     // BC value
	DragFunc DragFunction `json:"drag_func"` // FormFactor "G1" or "G7"
}

// Environment contains all the environmental variables
type Environment struct {
	Temperature        float64 `json:"temperature"`          // In Fahrenheit
	Altitude           int     `json:"altitude"`             // in ft
	Pressure           float64 `json:"pressure"`             // in Hg
	Humidity           float64 `json:"humidity"`             // Himidity percentage
	WindAngle          float64 `json:"wind_angle"`           // Wind direction in degrees (0 to 360)
	WindSpeed          float64 `json:"wind_speed"`           // Wind speed in mph
	PressureIsAbsolute bool    `json:"pressure_is_absolute"` // If true, only Pressure will be used
	Latitude           float64 `json:"latitude"`             // Latitude (deg)
	Azimuth            float64 `json:"azimuth"`              // Azimuth of fire (deg)
}

// Load includes data relevant to the load
type Load struct {
	Bullet         *Bullet `json:"bullet"`
	MuzzleVelocity float64 `json:"muzzle_velocity"` // in ft/s (default 3000)
}

// Rifle contains all the variables of the rifle
type Rifle struct {
	SightHeight        float64 `json:"sight_height"`         // Sight height in inches (default 1.5)
	BarrelTwist        float64 `json:"barrel_twist"`         // Barrel twist rate (default 7)
	TwistDirectionLeft bool    `json:"twist_direction_left"` // If twist direction is left
	ZeroRange          float64 `json:"zero_range"`           // Zero range of rifle in yrds (default 100)
}

// Solver is used to get solutions
type Solver struct {
	Load             *Load
	Rifle            *Rifle
	Environment      *Environment
	IncludeSpinDrift bool
	IncludeCoriolis  bool
}

// NewSolver creates a new Solver
func NewSolver(r *Rifle, e *Environment, l *Load, drift bool) *Solver {
	return &Solver{
		Rifle:            r,
		Environment:      e,
		Load:             l,
		IncludeSpinDrift: drift,
	}
}

// SetLoad sets the load
func (s *Solver) SetLoad(load *Load) { s.Load = load }

// SetEnvironment sets the environment
func (s *Solver) SetEnvironment(env *Environment) { s.Environment = env }

// SetRifle sets the rifle
func (s *Solver) SetRifle(rifle *Rifle) { s.Rifle = rifle }

// SetIncludeSpinDrift specifies where or not to include spin drift
func (s *Solver) SetIncludeSpinDrift(includeSpinDrift bool) { s.IncludeSpinDrift = includeSpinDrift }

// ErrMissingData - missing data error
var ErrMissingData = fmt.Errorf("missing rifle, environment, or target")

// ErrSolutionNotFound - no solution found error
var ErrSolutionNotFound = fmt.Errorf("no solution found")

// solve  is the basic solve function which takes a function
// telling when to stop (meaning solution found)
func (s *Solver) solve(stop func(*Solution) bool) (*Solution, error) {
	if s.Rifle == nil || s.Environment == nil || s.Load == nil {
		return nil, ErrMissingData
	}

	// Correct the DragCoefficient
	// for the Atmospheric conditions
	DragCoefficient := atmCorrect(
		s.Load.Bullet.BC.Value,
		float64(s.Environment.Altitude),
		s.Environment.Pressure,
		s.Environment.Temperature,
		s.Environment.Humidity,
	)
	if s.Environment.PressureIsAbsolute {
		DragCoefficient = pressureOnlyCorrect(
			s.Load.Bullet.BC.Value,
			s.Environment.Pressure,
		)
	}

	var (
		Vi              = s.Load.MuzzleVelocity
		WindSpeed       = s.Environment.WindSpeed
		WindAngle       = s.Environment.WindAngle
		DragFunc        = s.Load.Bullet.BC.DragFunc
		SightHeight     = s.Rifle.SightHeight
		ShootingAngle   = 0.0
		ZeroRange       = s.Rifle.ZeroRange
		ZeroAngle       = calcZeroAngle(DragFunc, DragCoefficient, Vi, SightHeight, ZeroRange, 0)
		StabilityFactor = gyroscopicStabilityFactor(s.Rifle, s.Load, s.Environment)
	)

	v := 0.0
	var vx, vx1, vy, vy1 float64
	var dv, dvx, dvy float64
	var x, y float64

	headwind := headWind(WindSpeed, WindAngle)
	crosswind := crossWind(WindSpeed, WindAngle)

	Gy := GRAVITY * math.Cos(degToRad((ShootingAngle + ZeroAngle)))
	Gx := GRAVITY * math.Sin(degToRad((ShootingAngle + ZeroAngle)))

	vx = Vi * math.Cos(degToRad(ZeroAngle))
	vy = Vi * math.Sin(degToRad(ZeroAngle))

	y = -SightHeight / 12

	t := 0.0
	for i := 0; i <= MAXRANGE; {

		vx1 = vx
		vy1 = vy
		v = math.Pow(math.Pow(vx, 2)+math.Pow(vy, 2), 0.5)
		dt := 0.5 / v

		// Compute acceleration using the drag function retardation
		dv = retard(DragFunc, DragCoefficient, v+headwind)
		dvx = -(vx / v) * dv
		dvy = -(vy / v) * dv

		// Compute velocity, including the resolved gravity vectors.
		vx = vx + dt*dvx + dt*Gx
		vy = vy + dt*dvy + dt*Gy

		wind := windage(crosswind, Vi, x, t+dt)
		drop := y * 12

		// Compute Spin drift and add to windage
		if s.IncludeSpinDrift {
			sd := spinDrift(t+dt, StabilityFactor)
			if s.Rifle.TwistDirectionLeft {
				sd = -sd
			}
			wind += sd
		}

		if s.IncludeCoriolis {
			wind += coriolisHorizontal(x/3, s.Environment.Latitude, t+dt)
			drop = coriolisVertical(
				s.Load.MuzzleVelocity,
				s.Environment.Latitude,
				s.Environment.Azimuth,
				drop,
			)
		}

		s := &Solution{
			Range:      i,
			RawRange:   x / 3,
			Drop:       drop,
			DropMOA:    -radToMOA(math.Atan(y / x)),
			Time:       t + dt,
			Windage:    wind,
			WindageMOA: radToMOA(math.Atan((wind / 12) / x)),
			Energy:     energy(s.Load.Bullet.Weight, v),
			Velocity:   v,
			VelocityX:  vx,
			VelocityY:  vy,
		}

		if stop(s) {
			return s, nil
		}

		if x/3 >= float64(i) {
			i++
		}

		// Compute position based on average velocity.
		x = x + dt*(vx+vx1)/2
		y = y + dt*(vy+vy1)/2

		if math.Abs(vy) > math.Abs(3*vx) {
			break
		}
		t = t + dt
	}

	return nil, ErrSolutionNotFound
}

// SolveSubsonic returns the solution where the
// bullet goes subsonic. You should include temperature
// in environment for this to be correct
func (s *Solver) SolveSubsonic() (*Solution, error) {
	temp := s.Environment.Temperature
	ss := 49.0223 * math.Pow(temp+459.67, 0.5)
	return s.solve(func(sol *Solution) bool {
		return sol.Velocity <= ss
	})
}

// SolveMaxRangeOGW uses the Matunas' Optimal Game Weight Formula to determine
// the maximum range for game of a given weight.
// It should be noted that this is an approximation
// - Takes the game weight in lbs
func (s *Solver) SolveMaxRangeOGW(weight float64) (*Solution, error) {
	return s.solve(func(sol *Solution) bool {
		ogw := math.Pow(sol.Velocity, 3.0) *
			math.Pow(s.Load.Bullet.Weight, 2.0) *
			1.5 * math.Pow(10.0, -12.0)
		return weight >= ogw
	})
}

// SolveFor uses the solver to solve for one specific range
func (s *Solver) SolveFor(rng int) (*Solution, error) {
	return s.solve(func(s *Solution) bool {
		return s.Range >= rng
	})
}

// Generate returns the solution set from range `min` to range `max`
// in increments of `inc`
func (s *Solver) Generate(min, max, inc int) []*Solution {
	// Let's not get crazy
	if max > MAXRANGE {
		max = MAXRANGE
	}
	if min < 0 {
		min = 0
	}
	if min > max {
		max, min = min, max
	}

	// Correct the DragCoefficient
	// for the Atmospheric conditions
	DragCoefficient := atmCorrect(
		s.Load.Bullet.BC.Value,
		float64(s.Environment.Altitude),
		s.Environment.Pressure,
		s.Environment.Temperature,
		s.Environment.Humidity,
	)
	if s.Environment.PressureIsAbsolute {
		DragCoefficient = pressureOnlyCorrect(
			s.Load.Bullet.BC.Value,
			s.Environment.Pressure,
		)
	}

	// Set up our variable
	var (
		Vi              = s.Load.MuzzleVelocity
		WindSpeed       = s.Environment.WindSpeed
		WindAngle       = s.Environment.WindAngle
		DragFunc        = s.Load.Bullet.BC.DragFunc
		SightHeight     = s.Rifle.SightHeight
		ShootingAngle   = 0.0 // TODO: remove, or add as an option
		ZeroRange       = s.Rifle.ZeroRange
		ZeroAngle       = calcZeroAngle(DragFunc, DragCoefficient, Vi, SightHeight, ZeroRange, 0)
		StabilityFactor = gyroscopicStabilityFactor(s.Rifle, s.Load, s.Environment)
	)

	// For now we will just append
	// we can maybe speed up by setting to (max-min)/inc
	Solutions := make([]*Solution, 0)

	// Set up integration variables
	v := 0.0
	var vx, vx1, vy, vy1 float64
	var dv, dvx, dvy float64
	var x, y float64

	headwind := headWind(WindSpeed, WindAngle)
	crosswind := crossWind(WindSpeed, WindAngle)

	// Gravitational accelerations
	Gy := GRAVITY * math.Cos(degToRad((ShootingAngle + ZeroAngle)))
	Gx := GRAVITY * math.Sin(degToRad((ShootingAngle + ZeroAngle)))

	// Break velocity into components
	vx = Vi * math.Cos(degToRad(ZeroAngle))
	vy = Vi * math.Sin(degToRad(ZeroAngle))

	// Bullet's starting height is
	// offset from the scope
	y = -SightHeight / 12

	t := 0.0
	for i := min; i <= max; {
		vx1 = vx
		vy1 = vy
		v = math.Pow(math.Pow(vx, 2)+math.Pow(vy, 2), 0.5)
		dt := 0.5 / v

		// Compute acceleration using the drag function retardation
		dv = retard(DragFunc, DragCoefficient, v+headwind)
		dvx = -(vx / v) * dv
		dvy = -(vy / v) * dv

		// Compute velocity, including the resolved gravity vectors.
		vx = vx + dt*dvx + dt*Gx
		vy = vy + dt*dvy + dt*Gy

		// Check if we've reached the next range
		// to add to the solution set
		if x/3 >= float64(i) {
			// Compute windage offset
			wind := windage(crosswind, Vi, x, t+dt)

			// Compute Spin drift and add to windage
			if s.IncludeSpinDrift {
				sd := spinDrift(t+dt, StabilityFactor)
				if s.Rifle.TwistDirectionLeft {
					sd = -sd
				}
				wind += sd
			}

			drop := y * 12
			if s.IncludeCoriolis {
				wind += coriolisHorizontal(x/3, s.Environment.Latitude, t+dt)
				drop = coriolisVertical(
					s.Load.MuzzleVelocity,
					s.Environment.Latitude,
					s.Environment.Azimuth,
					drop,
				)
			}

			// Calculations are done in ft
			// so some values need converted
			// for the solution
			Solutions = append(Solutions, &Solution{
				Range:      i,
				RawRange:   x / 3,                       // convert to yrds
				Drop:       drop,                        // convert to inches
				DropMOA:    -radToMOA(math.Atan(y / x)), // convert to MOA
				Time:       t + dt,
				Windage:    wind,
				WindageMOA: radToMOA(math.Atan((wind / 12) / x)), // convert to MOA
				Energy:     energy(s.Load.Bullet.Weight, v),
				Velocity:   v,
				VelocityX:  vx,
				VelocityY:  vy,
			})

			// Increment to next range
			i += inc
		}

		// Compute position based on average velocity.
		x = x + dt*(vx+vx1)/2
		y = y + dt*(vy+vy1)/2

		// If the bullet is dropping that fast
		// the computations are probably meaningless
		// (bullet is likely hitting terminal velocity)
		if math.Abs(vy) > math.Abs(3*vx) {
			break
		}
		t = t + dt
	}

	return Solutions
}

// SaveToCsv saves the solution set
// to a the csv file specified
func (s *Solver) SaveToCsv(path string, min, max, inc int) error {
	file, err := os.Create(path)
	if err != nil {
		return fmt.Errorf("could not create file: %v", err)
	}
	defer file.Close()

	csvWriter := csv.NewWriter(file)

	solutionSet := s.Generate(min, max, inc)

	// Add Header
	header := []string{
		"Range (yrds)",
		"Raw Range (yrds)",
		"Drop (in)",
		"Drop (MOA)",
		"Time (s)",
		"Wind (in)",
		"Wind (MOA)",
		"Energy (ft*lbs)",
		"Velocity (fps)",
		"VelocityX (fps)",
		"VelocityY (fps)",
	}
	csvWriter.Write(header)

	// Add Solution Set
	for _, s := range solutionSet {
		err := csvWriter.Write([]string{
			fmt.Sprintf("%d", s.Range),
			fmt.Sprintf("%0.2f", s.RawRange),
			fmt.Sprintf("%0.2f", s.Drop),
			fmt.Sprintf("%0.2f", s.DropMOA),
			fmt.Sprintf("%0.2f", s.Time),
			fmt.Sprintf("%0.2f", s.Windage),
			fmt.Sprintf("%0.2f", s.WindageMOA),
			fmt.Sprintf("%0.2f", s.Energy),
			fmt.Sprintf("%0.2f", s.Velocity),
			fmt.Sprintf("%0.2f", s.VelocityX),
			fmt.Sprintf("%0.2f", s.VelocityY),
		})
		if err != nil {
			return fmt.Errorf("could not write line: %v", err)
		}
	}

	csvWriter.Flush()

	if err := csvWriter.Error(); err != nil {
		return fmt.Errorf("could not flush csv: %v", err)
	}

	return nil
}

// ProjectedVelocity is used to calculate the projected muzzle velocity of a load
// if the bullet weight were changed, while it can be useful
// it is definitely just an estimate
func (l *Load) ProjectedVelocity(weight float64) float64 {
	return l.MuzzleVelocity * math.Pow(l.Bullet.Weight/weight, 2.0)
}

/* Functions for correcting for atmosphere. */

// Correct for Temperature + humidity
func calcFR(Temperature, Pressure, RelativeHumidity float64) float64 {
	VPw := 4e-6*math.Pow(Temperature, 3) - 0.0004*math.Pow(Temperature, 2) + 0.0234*Temperature - 0.2517
	FRH := 0.995 * (Pressure / (Pressure - (0.3783)*(RelativeHumidity)*VPw))
	return FRH
}

// Correct for different Pressure than standard
func calcFP(Pressure float64) float64 {
	Pstd := 29.53 // in-hg
	return (Pressure - Pstd) / (Pstd)
}

// Correct at Temperature + Altitude
func calcFT(Temperature, Altitude float64) float64 {
	Tstd := -0.0036*Altitude + 59
	FT := (Temperature - Tstd) / (459.6 + Tstd)
	return FT
}

// Correct for Altitude
func calcFA(Altitude float64) float64 {
	fa := -4e-15*math.Pow(Altitude, 3) + 4e-10*math.Pow(Altitude, 2) - 3e-5*Altitude + 1
	return (1 / fa)
}

// Combine all the atmospheric corrections
// and return the corrected BC
func atmCorrect(dragCoefficient, altitude, pressure, temperature, relativeHumidity float64) float64 {
	FA := calcFA(altitude)
	FT := calcFT(temperature, altitude)
	FR := calcFR(temperature, pressure, relativeHumidity)
	FP := calcFP(pressure)

	// Calculate the atmospheric correction factor
	CD := (FA * (1 + FT - FP) * FR)
	return dragCoefficient * CD
}

// Only correct for pressure difference (when pressure is absolute)
func pressureOnlyCorrect(dragCoefficient, pressure float64) float64 {
	FP := calcFP(pressure)
	return dragCoefficient * (1 - FP)
}

// Amount of horizontal coriolis correction
// Needs range, latitude and time of flight
func coriolisHorizontal(rng, lat, t float64) float64 {
	deflection := EARTH_ROTATION * rng * math.Sin(degToRad(lat)) * t
	// Convert to inches and return
	return deflection * 12
}

// Vertical correctionf or coriolis
// Needs MuzzleVelocity, Latitude, Azimuth, Drop
func coriolisVertical(mv, lat, az, drop float64) float64 {
	// Calculate the correction factor
	cf := 1 - ((2*EARTH_ROTATION*mv)/GRAVITY)*math.Cos(degToRad(lat))*math.Sin(degToRad(az))

	return drop * cf
}

// Calculate Spin Drift
// Takes time of flight (s) and gyroscopic stability factor
// Returns amount of drift (in)
func spinDrift(t, sf float64) float64 {
	return 1.25 * (sf + 1.2) * math.Pow(t, 1.83)
}

// Calculates the gyroscopic stability factor
// Using the Miller Stability Formula and
// Velocity correction for high muzzle velocities
func gyroscopicStabilityFactor(r *Rifle, l *Load, e *Environment) float64 {
	sf := millerStabilityFormula(l.Bullet.Weight, r.BarrelTwist, l.Bullet.Caliber, l.Bullet.Length)
	if l.MuzzleVelocity > 2800 {
		sf = millerStabilityFormulaVelocityCorrection(sf, l.MuzzleVelocity)
	}
	sf = millerStabilityFormulaAtmosphericCorrection(sf, e.Temperature, e.Pressure)
	return sf
}

// MillerStabilityFormula takes
// m: mass of bullet (gr)
// t: rifling twist rate
// d: diameter/caliber of bullet
// l: length of bullet (in)
// and returns the gyroscopic stability factor
func millerStabilityFormula(m, t, d, l float64) float64 {
	return (30.0 * m) / (math.Pow(t/d, 2.0) * math.Pow(d, 3.0) * (l / d) * (1.0 + math.Pow(l/d, 2.0)))
}

// The basic MillerStabilityFormual is valid for muzzle velocities
// up to 2800fps. For higher muzzle velocities
// This function corrects for higher muzzle velocities
func millerStabilityFormulaVelocityCorrection(sf, mv float64) float64 {
	return sf * math.Pow(mv/2800, 1.0/3.0)
}

// Atmospheric correction factor for MillerStability Formula
func millerStabilityFormulaAtmosphericCorrection(sf, temp, pressure float64) float64 {
	return sf * (((temp + 460) / 519.0) * (29.92 / pressure))
}

// Calculate Energy
// takes mass(gr) and velocity(ft/s)
func energy(mass, velocity float64) float64 {
	return (mass * velocity * velocity) / 450800
}

// Calculate windage
// Wind Deflection = 17.6 * CrosswindSpeed * LagTime
// LagTime is (actual time of flight) - (time of flight in a vacuum)
func windage(WindSpeed, mv, rng, t float64) float64 {
	Vw := WindSpeed * 17.60
	return Vw * (t - rng/mv)
}

// Headwind is positive at WindAngle=0
func headWind(WindSpeed, WindAngle float64) float64 {
	Wangle := degToRad(WindAngle)
	return (math.Cos(Wangle) * WindSpeed)
}

// Positive is from Shooter's Right to Left (Wind from 90 degree)
func crossWind(WindSpeed, WindAngle float64) float64 {
	Wangle := degToRad(WindAngle)
	return (math.Sin(Wangle) * WindSpeed)
}

// CalcZeroAngle computes the angle
// that would have been used to zero the
// rifle at the given range
func calcZeroAngle(DragFunc DragFunction, DragCoefficient, Vi, SightHeight, ZeroRange, yIntercept float64) float64 {

	// Numerical Integration variables
	dt := 1 / Vi       // The solution accuracy generally doesn't suffer if its within a foot for each second of time.
	var dAngle float64 // The change in the bore angle used to iterate in on the correct zero angle.

	// State variables for each integration loop.
	var t, x float64
	var y = -SightHeight / 12.0
	var v, vx, vy float64    // velocity
	var vx1, vy1 float64     // Last frame's velocity, used for computing average velocity.
	var dv, dvx, dvy float64 // acceleration
	var Gx, Gy float64       // Gravitational acceleration

	var angle float64 // The actual angle of the bore.

	var quit bool // We know it's time to quit our successive approximation loop when this is true.

	// Start with a very coarse angular change, to quickly solve even large launch angle problems.
	dAngle = degToRad(14.0)

	// The general idea here is to start at 0 degrees elevation, and increase the elevation by 14 degrees
	// until we are above the correct elevation.  Then reduce the angular change by half, and begin reducing
	// the angle.  Once we are again below the correct angle, reduce the angular change by half again, and go
	// back up.  This allows for a fast successive approximation of the correct elevation, usually within less
	// than 20 iterations.
	for angle = 0.0; !quit; angle = angle + dAngle {
		vy = Vi * math.Sin(angle)
		vx = Vi * math.Cos(angle)
		Gx = GRAVITY * math.Sin(angle)
		Gy = GRAVITY * math.Cos(angle)

		for t, x, y = 0.0, 0.0, -SightHeight/12; x <= ZeroRange*3; t = t + dt {
			vy1 = vy
			vx1 = vx
			v = math.Pow((math.Pow(vx, 2) + math.Pow(vy, 2)), 0.5)
			dt = 1 / v

			dv = retard(DragFunc, DragCoefficient, v)
			dvy = -dv * vy / v * dt
			dvx = -dv * vx / v * dt

			vx = vx + dvx
			vy = vy + dvy
			vy = vy + dt*Gy
			vx = vx + dt*Gx

			x = x + dt*(vx+vx1)/2
			y = y + dt*(vy+vy1)/2
			// Break early to save CPU time if we won't find a solution.
			if vy < 0 && y < yIntercept {
				break
			}
			if vy > 3*vx {
				break
			}
		}

		if y > yIntercept && dAngle > 0 {
			dAngle = -dAngle / 2
		}
		if y < yIntercept && dAngle < 0 {
			dAngle = -dAngle / 2
		}

		// If our accuracy is sufficient, we can stop approximating.
		if math.Abs(dAngle) < moaToRad(0.01) {
			quit = true
		}
		// If we exceed the 45 degree launch angle
		// then the projectile just won't get there
		// so we stop trying.
		if angle > degToRad(45) {
			quit = true
		}

	}

	return radToDeg(angle) // Convert to degrees for return value.
}

/* Angular conversion functions */
func degToMOA(deg float64) float64 {
	return deg * 60
}
func degToRad(deg float64) float64 {
	return deg * math.Pi / 180
}
func moaToDeg(moa float64) float64 {
	return moa / 60
}
func moaToRad(moa float64) float64 {
	return moa / 60 * math.Pi / 180
}
func radToDeg(rad float64) float64 {
	return rad * 180 / math.Pi
}
func radToMOA(rad float64) float64 {
	return rad * 60 * 180 / math.Pi
}

// Retard corrects the velocity for ballistic drag.
func retard(DragFunc DragFunction, DragCoefficient, Velocity float64) float64 {
	vp := Velocity
	A := -1.0
	M := -1.0

	switch DragFunc {
	case G1:
		if vp > 4230 {
			A = 1.477404177730177e-04
			M = 1.9565
		} else if vp > 3680 {
			A = 1.920339268755614e-04
			M = 1.925
		} else if vp > 3450 {
			A = 2.894751026819746e-04
			M = 1.875
		} else if vp > 3295 {
			A = 4.349905111115636e-04
			M = 1.825
		} else if vp > 3130 {
			A = 6.520421871892662e-04
			M = 1.775
		} else if vp > 2960 {
			A = 9.748073694078696e-04
			M = 1.725
		} else if vp > 2830 {
			A = 1.453721560187286e-03
			M = 1.675
		} else if vp > 2680 {
			A = 2.162887202930376e-03
			M = 1.625
		} else if vp > 2460 {
			A = 3.209559783129881e-03
			M = 1.575
		} else if vp > 2225 {
			A = 3.904368218691249e-03
			M = 1.55
		} else if vp > 2015 {
			A = 3.222942271262336e-03
			M = 1.575
		} else if vp > 1890 {
			A = 2.203329542297809e-03
			M = 1.625
		} else if vp > 1810 {
			A = 1.511001028891904e-03
			M = 1.675
		} else if vp > 1730 {
			A = 8.609957592468259e-04
			M = 1.75
		} else if vp > 1595 {
			A = 4.086146797305117e-04
			M = 1.85
		} else if vp > 1520 {
			A = 1.954473210037398e-04
			M = 1.95
		} else if vp > 1420 {
			A = 5.431896266462351e-05
			M = 2.125
		} else if vp > 1360 {
			A = 8.847742581674416e-06
			M = 2.375
		} else if vp > 1315 {
			A = 1.456922328720298e-06
			M = 2.625
		} else if vp > 1280 {
			A = 2.419485191895565e-07
			M = 2.875
		} else if vp > 1220 {
			A = 1.657956321067612e-08
			M = 3.25
		} else if vp > 1185 {
			A = 4.745469537157371e-10
			M = 3.75
		} else if vp > 1150 {
			A = 1.379746590025088e-11
			M = 4.25
		} else if vp > 1100 {
			A = 4.070157961147882e-13
			M = 4.75
		} else if vp > 1060 {
			A = 2.938236954847331e-14
			M = 5.125
		} else if vp > 1025 {
			A = 1.228597370774746e-14
			M = 5.25
		} else if vp > 980 {
			A = 2.916938264100495e-14
			M = 5.125
		} else if vp > 945 {
			A = 3.855099424807451e-13
			M = 4.75
		} else if vp > 905 {
			A = 1.185097045689854e-11
			M = 4.25
		} else if vp > 860 {
			A = 3.566129470974951e-10
			M = 3.75
		} else if vp > 810 {
			A = 1.045513263966272e-08
			M = 3.25
		} else if vp > 780 {
			A = 1.291159200846216e-07
			M = 2.875
		} else if vp > 750 {
			A = 6.824429329105383e-07
			M = 2.625
		} else if vp > 700 {
			A = 3.569169672385163e-06
			M = 2.375
		} else if vp > 640 {
			A = 1.839015095899579e-05
			M = 2.125
		} else if vp > 600 {
			A = 5.71117468873424e-05
			M = 1.950
		} else if vp > 550 {
			A = 9.226557091973427e-05
			M = 1.875
		} else if vp > 250 {
			A = 9.337991957131389e-05
			M = 1.875
		} else if vp > 100 {
			A = 7.225247327590413e-05
			M = 1.925
		} else if vp > 65 {
			A = 5.792684957074546e-05
			M = 1.975
		} else if vp > 0 {
			A = 5.206214107320588e-05
			M = 2.000
		}
		break
	case G2:
		if vp > 1674 {
			A = 0.0079470052136733
			M = 1.36999902851493
		} else if vp > 1172 {
			A = 1.00419763721974e-03
			M = 1.65392237010294
		} else if vp > 1060 {
			A = 7.15571228255369e-23
			M = 7.91913562392361
		} else if vp > 949 {
			A = 1.39589807205091e-10
			M = 3.81439537623717
		} else if vp > 670 {
			A = 2.34364342818625e-04
			M = 1.71869536324748
		} else if vp > 335 {
			A = 1.77962438921838e-04
			M = 1.76877550388679
		} else if vp > 0 {
			A = 5.18033561289704e-05
			M = 1.98160270524632
		}
		break
	case G5:
		if vp > 1730 {
			A = 7.24854775171929e-03
			M = 1.41538574492812
		} else if vp > 1228 {
			A = 3.50563361516117e-05
			M = 2.13077307854948
		} else if vp > 1116 {
			A = 1.84029481181151e-13
			M = 4.81927320350395
		} else if vp > 1004 {
			A = 1.34713064017409e-22
			M = 7.8100555281422
		} else if vp > 837 {
			A = 1.03965974081168e-07
			M = 2.84204791809926
		} else if vp > 335 {
			A = 1.09301593869823e-04
			M = 1.81096361579504
		} else if vp > 0 {
			A = 3.51963178524273e-05
			M = 2.00477856801111
		}
		break
	case G6:
		if vp > 3236 {
			A = 0.0455384883480781
			M = 1.15997674041274
		} else if vp > 2065 {
			A = 7.167261849653769e-02
			M = 1.10704436538885
		} else if vp > 1311 {
			A = 1.66676386084348e-03
			M = 1.60085100195952
		} else if vp > 1144 {
			A = 1.01482730119215e-07
			M = 2.9569674731838
		} else if vp > 1004 {
			A = 4.31542773103552e-18
			M = 6.34106317069757
		} else if vp > 670 {
			A = 2.04835650496866e-05
			M = 2.11688446325998
		} else if vp > 0 {
			A = 7.50912466084823e-05
			M = 1.92031057847052
		}
		break
	case G7:
		if vp > 4200 {
			A = 1.29081656775919e-09
			M = 3.24121295355962
		} else if vp > 3000 {
			A = 0.0171422231434847
			M = 1.27907168025204
		} else if vp > 1470 {
			A = 2.33355948302505e-03
			M = 1.52693913274526
		} else if vp > 1260 {
			A = 7.97592111627665e-04
			M = 1.67688974440324
		} else if vp > 1110 {
			A = 5.71086414289273e-12
			M = 4.3212826264889
		} else if vp > 960 {
			A = 3.02865108244904e-17
			M = 5.99074203776707
		} else if vp > 670 {
			A = 7.52285155782535e-06
			M = 2.1738019851075
		} else if vp > 540 {
			A = 1.31766281225189e-05
			M = 2.08774690257991
		} else if vp > 0 {
			A = 1.34504843776525e-05
			M = 2.08702306738884
		}
		break
	case G8:
		if vp > 3571 {
			A = .0112263766252305
			M = 1.33207346655961
		} else if vp > 1841 {
			A = .0167252613732636
			M = 1.28662041261785
		} else if vp > 1120 {
			A = 2.20172456619625e-03
			M = 1.55636358091189
		} else if vp > 1088 {
			A = 2.0538037167098e-16
			M = 5.80410776994789
		} else if vp > 976 {
			A = 5.92182174254121e-12
			M = 4.29275576134191
		} else if vp > 0 {
			A = 4.3917343795117e-05
			M = 1.99978116283334
		}
		break
	default:
		break

	}

	if A != -1 && M != -1 && vp > 0 && vp < 10000 {
		val := A * math.Pow(vp, M) / DragCoefficient
		return val
	}
	return -1
}
