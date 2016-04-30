using UnityEngine;

namespace Barracuda
{
    public delegate float EasingMode(float x, float t, float b, float c, float d);

    public static class Easing
    {
        /// <summary>
        /// Slerp the specified easingMode, time and duration.
        /// if time is not in the range of 0 to 1, the value will be rounded into 0 or 1.
        /// </summary>
        public static float Slerp(EasingMode easingMode, float time, float duration) {
            if (time < 0) {
                return 0;
            }
            if (time > duration) {
                return 1;
            }
            return easingMode(time, time, 0, 1, duration);
        }

        public static EasingMode FromAnimationCurve(AnimationCurve curve)
        {
            return new EasingMode(delegate(float x, float t, float b, float c, float d)
            {
                return curve.Evaluate(x / d);
            });
        }

        public static EasingMode EaseOut
        {
            get
            {
                return Easing.EaseOutQuad;
            }
        }

        public static EasingMode EaseIn
        {
            get
            {
                return Easing.EaseInQuad;
            }
        }

        public static float Linear(float x, float t, float b, float c, float d)
        {
            return t / d;
        }

        public static float EaseInQuad(float x, float t, float b, float c, float d)
        {
            return c * (t /= d) * t + b;
        }

        public static float EaseOutQuad(float x, float t, float b, float c, float d)
        {
            return -c * (t /= d) * (t - 2) + b;
        }

        public static float EaseInOutQuad(float x, float t, float b, float c, float d)
        {
            if ((t /= d / 2) < 1)
            {
                return c / 2 * t * t + b;
            }
            return -c / 2 * ((--t) * (t - 2) - 1) + b;
        }

        public static float EaseInCubic(float x, float t, float b, float c, float d)
        {
            return c * (t /= d) * t * t + b;
        }

        public static float EaseOutCubic(float x, float t, float b, float c, float d)
        {
            return c * ((t = t / d - 1) * t * t + 1) + b;
        }

        public static float EaseInOutCubic(float x, float t, float b, float c, float d)
        {
            if ((t /= d / 2) < 1)
            {
                return c / 2 * t * t * t + b;
            }
            return c / 2 * ((t -= 2) * t * t + 2) + b;
        }

        public static float EaseInQuart(float x, float t, float b, float c, float d)
        {
            return c*(t/=d)*t*t*t + b;
        }

        public static float EaseOutQuart(float x, float t, float b, float c, float d)
        {
            return -c * ((t=t/d-1)*t*t*t - 1) + b;
        }

        public static float EaseInOutQuart(float x, float t, float b, float c, float d)
        {
            if ((t/=d/2) < 1) return c/2*t*t*t*t + b;
            return -c/2 * ((t-=2)*t*t*t - 2) + b;
        }

        public static float EaseInQuint(float x, float t, float b, float c, float d)
        {
            return c*(t/=d)*t*t*t*t + b;
        }

        public static float EaseOutQuint(float x, float t, float b, float c, float d)
        {
            return c*((t=t/d-1)*t*t*t*t + 1) + b;
        }

        public static float EaseInOutQuint(float x, float t, float b, float c, float d)
        {
            if ((t/=d/2) < 1) return c/2*t*t*t*t*t + b;
            return c/2*((t-=2)*t*t*t*t + 2) + b;
        }

        public static float EaseInSine(float x, float t, float b, float c, float d)
        {
            return -c * Mathf.Cos(t/d * (Mathf.PI/2)) + c + b;
        }

        public static float EaseOutSine(float x, float t, float b, float c, float d)
        {
            return c * Mathf.Sin(t/d * (Mathf.PI/2)) + b;
        }

        public static float EaseInOutSine(float x, float t, float b, float c, float d)
        {
            return -c/2 * (Mathf.Cos(Mathf.PI*t/d) - 1) + b;
        }

        public static float EaseInExpo(float x, float t, float b, float c, float d)
        {
            return (t==0) ? b : c * Mathf.Pow(2, 10 * (t/d - 1)) + b;
        }

        public static float EaseOutExpo(float x, float t, float b, float c, float d)
        {
            return (t==d) ? b+c : c * (-Mathf.Pow(2, -10 * t/d) + 1) + b;
        }

        public static float EaseInOutExpo(float x, float t, float b, float c, float d)
        {
            if (t==0) return b;
            if (t==d) return b+c;
            if ((t/=d/2) < 1) return c/2 * Mathf.Pow(2, 10 * (t - 1)) + b;
            return c/2 * (-Mathf.Pow(2, -10 * --t) + 2) + b;
        }

        public static float EaseInCirc(float x, float t, float b, float c, float d)
        {
            return -c * (Mathf.Sqrt(1 - (t/=d)*t) - 1) + b;
        }

        public static float EaseOutCirc(float x, float t, float b, float c, float d)
        {
            return c * Mathf.Sqrt(1 - (t=t/d-1)*t) + b;
        }

        public static float EaseInOutCirc(float x, float t, float b, float c, float d)
        {
            if ((t/=d/2) < 1) return -c/2 * (Mathf.Sqrt(1 - t*t) - 1) + b;
            return c/2 * (Mathf.Sqrt(1 - (t-=2)*t) + 1) + b;
        }

        public static float EaseInElastic(float x, float t, float b, float c, float d)
        {
            var s = 1.70158f;
            var p = d * .3f;
            var a = c;
            if (t == 0) return b;
            if ((t/=d) == 1) return b+c;
            if (a < Mathf.Abs(c))
            {
                a = c;
                s = p/4;
            }
            else
            {
                s = p / (2 * Mathf.PI) * Mathf.Asin(c / a);
            }
            return -(a * Mathf.Pow(2,10 * (t -= 1)) * Mathf.Sin((t * d - s)*(2 * Mathf.PI) / p)) + b;
        }

        public static float EaseOutElastic(float x, float t, float b, float c, float d)
        {
            var s = 1.70158f;
            var p = d * .3f;
            var a = c;
            if (t == 0) return b;
            if ((t /= d) == 1) return b+c;
            if (a < Mathf.Abs(c)) {
                a = c;
                s = p / 4;
            }
            else
            {
                s = p / (2 * Mathf.PI) * Mathf.Asin(c / a);
            }
            return a * Mathf.Pow(2,-10 * t) * Mathf.Sin((t * d - s) * (2 * Mathf.PI) / p) + c + b;
        }


    	public static float VibratePendulum(float x, float t, float b, float c, float d)
    	{
    		var p = (t / d);
    		return Mathf.Sin (4.0f * p * Mathf.PI) * (1 - p);
    	}

        public static float EaseInOutElastic(float x, float t, float b, float c, float d)
        {
            var s = 1.70158f;
            var p = 0f;
            var a = c;
            if (t == 0)
            {
                return b;
            }
            if ((t /= d / 2) == 2)
            {
                return b + c;
            }
            if (p != 0)
            {
                p = d * (.3f * 1.5f);
            }
            if (a < Mathf.Abs(c)) {
                a = c;
                s = p / 4;
            }
            else
            {
                s = p / (2 * Mathf.PI) * Mathf.Asin(c / a);
            }
            if (t < 1)
            {
                return -.5f * (a * Mathf.Pow(2,10 * (t -= 1)) * Mathf.Sin((t * d - s) * (2 * Mathf.PI) / p)) + b;
            }
            return a * Mathf.Pow(2,-10 * (t -= 1)) * Mathf.Sin((t * d - s) * (2 * Mathf.PI) / p ) * .5f + c + b;
        }

        public static float EaseInBack(float x, float t, float b, float c, float d)
        {
            var s = 1.70158f;
            return c * (t /= d) * t * ((s + 1) * t - s) + b;
        }

        public static float EaseOutBack(float x, float t, float b, float c, float d)
        {
            var s = 1.70158f;
            return c * ((t = (t/d - 1)) * t * ((s + 1) * t + s) + 1) + b;
        }

        public static float EaseInOutBack(float x, float t, float b, float c, float d)
        {
            var s = 1.70158f;
            if ((t /= (d / 2)) < 1)
            {
                return c / 2 * (t * t * (((s *= (1.525f)) + 1) * t - s)) + b;
            }
            return c / 2 * ((t -= 2) * t * (((s *= (1.525f)) + 1) * t + s) + 2) + b;
        }

        public static float EaseInBounce(float x, float t, float b, float c, float d)
        {
            return c - Easing.EaseOutBounce (x, d - t, 0, c, d) + b;
        }

        public static float EaseOutBounce(float x, float t, float b, float c, float d)
        {
            if ((t /= d) < (1 / 2.75f))
            {
                return c * (7.5625f * t * t) + b;
            }
            else if (t < (2/2.75f))
            {
                return c * (7.5625f * (t -= (1.5f/2.75f)) * t + .75f) + b;
            }
            else if (t < (2.5f / 2.75f))
            {
                return c * (7.5625f * (t -= (2.25f / 2.75f)) * t + .9375f) + b;
            }
            else
            {
                return c * (7.5625f * (t -= (2.625f / 2.75f)) * t + .984375f) + b;
            }
        }

        public static float EaseInOutBounce(float x, float t, float b, float c, float d)
        {
            if (t < (d / 2))
            {
                return Easing.EaseInBounce (x, t * 2, 0, c, d) * .5f + b;
            }
            return Easing.EaseOutBounce (x, (t * 2) - d, 0, c, d) * .5f + (c * .5f) + b;
        }

    }
}
