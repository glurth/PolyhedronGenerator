using Cysharp.Threading.Tasks;

namespace EyE.Threading
{
    /// <summary>
    /// provide a way to reference a single/the same bool- stored in here- regardless of thread context
    /// </summary>
    public class CancelBoolRef
    {
        public volatile bool doCancel = false;
    }

    /// <summary>
    /// Utility class to help throttle async task execution by yielding control
    /// back to Unity after a specified time slice has elapsed.
    /// Set the runSynchronously constructor param to override, and never actually yield.  (this allows calling code to stay the same, but handle both sync and async)
    /// </summary>
    public class YieldTimer
    {
        private readonly int timeSlice; // Time slice in milliseconds
        private readonly System.Diagnostics.Stopwatch timer = new System.Diagnostics.Stopwatch();
        private readonly CancelBoolRef cancelRef; // Optional cancellation reference
        private readonly bool runSynchronously;

        public YieldTimer(int timeSlice = 10, CancelBoolRef cancelRef = null, bool runSynchronously = false)
        {
            this.timeSlice = timeSlice;
            this.cancelRef = cancelRef;
            this.runSynchronously = runSynchronously;
            if (!runSynchronously) timer.Start();
        }
        public YieldTimer(CancelBoolRef cancelRef, bool runSynchronously = false)
        {
            this.timeSlice = 10;
            this.cancelRef = cancelRef;
            this.runSynchronously = runSynchronously;
            if (!runSynchronously) timer.Start();
        }
        /// <summary>
        /// If the processing timer has reached or exceeded the allotted time-slice, will yield until the next unity update.
        /// </summary>
        /// <returns></returns>
        public async UniTask YieldOnTimeSlice()
        {
            if (runSynchronously) return;
            if (cancelRef != null && cancelRef.doCancel)
                throw new System.OperationCanceledException();
            if (timer.ElapsedMilliseconds > timeSlice)
            {
                await UniTask.Yield(PlayerLoopTiming.Update);
                timer.Restart(); // Restart timer after yielding
            }
        }
    }
}